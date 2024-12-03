from collections import defaultdict
import numpy as np
import pandas as pd
import torch
import torch.nn.functional as F
import scipy
import scipy.io
from sklearn.preprocessing import label_binarize
import torch_geometric.transforms as T

# 因为修改了data_utils文件，直接import得到的一直是之前未修改的代码，所以要reload一下，先import，再reload
import data_utils
import importlib
importlib.reload(data_utils)

from data_utils import rand_train_test_idx, even_quantile_labels, to_sparse_tensor, dataset_drive_url, class_rand_splits, rand_divide_idx

from torch_geometric.datasets import Planetoid, Amazon, Coauthor
from torch_geometric.utils import degree
import os

from google_drive_downloader import GoogleDriveDownloader as gdd

import networkx as nx
import scipy.sparse as sp

from ogb.nodeproppred import NodePropPredDataset

import scanpy as sc
from sklearn.feature_extraction.text import TfidfTransformer


class NCDataset(object):
    def __init__(self, name):
        """
        based off of ogb NodePropPredDataset
        https://github.com/snap-stanford/ogb/blob/master/ogb/nodeproppred/dataset.py
        Gives torch tensors instead of numpy arrays
            - name (str): name of the dataset
            - root (str): root directory to store the dataset folder
            - meta_dict: dictionary that stores all the meta-information about data. Default is None, 
                    but when something is passed, it uses its information. Useful for debugging for external contributers.

        Usage after construction:

        split_idx = dataset.get_idx_split()
        train_idx, valid_idx, test_idx = split_idx["train"], split_idx["valid"], split_idx["test"]
        graph, label = dataset[0]

        Where the graph is a dictionary of the following form:
        dataset.graph = {'edge_index': edge_index,
                         'edge_feat': None,
                         'node_feat': node_feat,
                         'num_nodes': num_nodes}
        For additional documentation, see OGB Library-Agnostic Loader https://ogb.stanford.edu/docs/nodeprop/

        """

        self.name = name  # original name, e.g., ogbn-proteins
        self.graph = {}
        self.label = None
        self.metadata = [] # 用于记录数据来源，方便后续求loss或是数据集划分
        self.n_infer = 0 # 用于记录edge中锚点对的个数
        self.n_infer_intra1 = 0
        self.n_data1 = 0 # index从0开始计数，前面的是data1，后面的是data2
        self.n_data2 = 0
        self.num_celltype = {}
        self.edge = None # 用于记载边

    def get_idx_split(self, split_type='random', train_prop=.8, valid_prop=.1, label_num_per_class=20):
        """
        划分数据集，返回索引
        split_type: 'random' for random splitting, 'class' for splitting with equal node num per class
        train_prop: The proportion of dataset for train split. Between 0 and 1.
        这里我就只先改出了一个分割模式，函数加到了data_utils.py里面
        
        valid_prop: The proportion of dataset for validation split. Between 0 and 1.
        label_num_per_class: num of nodes per class
        """

        if split_type == 'random':
            train_idx, valid_idx, test_idx, query_idx = rand_divide_idx(
                self.label, self.metadata, train_prop=train_prop, valid_prop=valid_prop)
            split_idx = {'train': train_idx,
                         'valid': valid_idx,
                         'test': test_idx,
                         'query': query_idx}
        elif split_type == 'class':
            print('未定义！')
        return split_idx

    def __getitem__(self, idx):
        assert idx == 0, 'This dataset has only one graph'
        return self.graph, self.label

    def __len__(self):
        return 1

    def __repr__(self):
        return '{}({})'.format(self.__class__.__name__, len(self))

def input_dataset(data_dir, name):
    # data_dir默认为'../data/'，name为自定义的数据集名称
    dataset = NCDataset('dataset')

    DataPath1 = '{}/rna.h5ad'.format(data_dir)
    DataPath2 = '{}/atac.h5ad'.format(data_dir)
    LabelsPath1 = '{}label1.csv'.format(data_dir)
    LabelsPath2 = '{}label2.csv'.format(data_dir)
    
    inter_graphPath = '{}inter_graph.csv'.format(data_dir)
    intra_graphPath1 = '{}intra_graph1.csv'.format(data_dir)
    intra_graphPath2 = '{}intra_graph2.csv'.format(data_dir)

    #' read the data
    data1 = sc.read_h5ad(DataPath1)
    data2 = sc.read_h5ad(DataPath2)
    lab_label1 = pd.read_csv(LabelsPath1, header=0, index_col=False, sep=',')
    lab_label2 = pd.read_csv(LabelsPath2, header=0, index_col=False, sep=',')
    inter_graph = pd.read_csv(inter_graphPath, index_col=0, sep=',').to_numpy()
    intra_graph1 = pd.read_csv(intra_graphPath1, index_col=0, sep=',').to_numpy()
    intra_graph2 = pd.read_csv(intra_graphPath2, index_col=0, sep=',').to_numpy()
    
    if isinstance(data1.X, sp.csr_matrix):
        data1.X = data1.X.toarray()
    if isinstance(data2.X, sp.csr_matrix):
        data2.X = data2.X.toarray()
    
    adata = data1
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.scale(adata)
    data1 = adata.X

    adata = data2
    # sc.pp.normalize_total(adata, target_sum=1e4)
    # sc.pp.log1p(adata)
    tfidf = TfidfTransformer()
    adata.X = tfidf.fit_transform(adata.X).toarray()
    sc.pp.scale(adata)
    data2 = adata.X

    data1 = torch.from_numpy(data1)
    data2 = torch.from_numpy(data2)

    features = torch.cat((data1, data2), dim=0).to(torch.float32)

    lab_label1 = pd.read_csv(LabelsPath1, header=0, index_col=False, sep=',')
    lab_label2 = pd.read_csv(LabelsPath2, header=0, index_col=False, sep=',')
    lab_label1.columns = ['type']
    lab_label2.columns = ['type']
    label = pd.concat([lab_label1, lab_label2], ignore_index=True)
    types = np.unique(label['type']).tolist()
    label = label.values.tolist()
    label = [item[0] for item in label] # 把二维list转化为一维list
    num_celltype={} # 记录的celltype与数字的对应关系
    for i,j in enumerate(types):
        label = [i if x == j else x for x in label]
        num_celltype.update({i : j})
    

    inter_graph = torch.from_numpy(inter_graph)
    intra_graph1 = torch.from_numpy(intra_graph1)
    intra_graph2 = torch.from_numpy(intra_graph2)
    edge_index = torch.cat((inter_graph, intra_graph1, intra_graph2), dim=0) # 
    edge_index = torch.transpose(edge_index, 0, 1)
    dataset.edge = torch.cat((inter_graph, intra_graph2), dim=0)
    dataset.edge = torch.transpose(dataset.edge, 0, 1)

    metadata = [1] * len(lab_label1) + [2] * len(lab_label2)
    
    dataset.n_infer = inter_graph.shape[0]
    dataset.n_infer_intra1 = torch.cat((inter_graph, intra_graph1),dim=0).shape[0]
    
    dataset.n_data1 = data1.shape[0]
    dataset.n_data2 = data2.shape[0]

    dataset.graph = {'edge_index': edge_index,
                     'edge_feat': None,
                     'node_feat': features,
                     'num_nodes': len(label)}
    dataset.label = torch.LongTensor(label)
    dataset.metadata = torch.LongTensor(metadata)
    dataset.num_celltype = num_celltype
    
    return dataset

