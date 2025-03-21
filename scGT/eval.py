import torch
import torch.nn.functional as F
import torch.nn as nn
import numpy as np
from sklearn.metrics import silhouette_score,silhouette_samples, accuracy_score, f1_score
import pandas as pd
from sklearn.neighbors import NearestNeighbors

def eval_acc(y_true, y_pred):
    acc_list = []
    y_true = y_true.detach().cpu().numpy()
    y_pred = y_pred.argmax(dim=-1, keepdim=True).detach().cpu().numpy()

    for i in range(y_true.shape[1]):
        is_labeled = y_true[:, i] == y_true[:, i]
        correct = y_true[is_labeled, i] == y_pred[is_labeled, i]
        acc_list.append(float(np.sum(correct)) / len(correct))
    
    return sum(acc_list) / len(acc_list)

@torch.no_grad()
def evaluate_cpu(model, dataset, split_idx, eval_func, criterion, args, result=None):
    model.eval()

    model.to(torch.device("cpu"))
    dataset.label = dataset.label.to(torch.device("cpu"))
    adjs_, x = dataset.graph['adjs'], dataset.graph['node_feat']
    adjs = []
    adjs.append(adjs_[0])
    for k in range(1):
        adjs.append(adjs_[k + 1])
    out, _, z1, z2 = model(x, adjs)

    train_acc = eval_func(
        dataset.label[split_idx['train']], out[split_idx['train']])

    return train_acc

def silhouette(embedding_umap, label, tech):
    ## 轮廓系数
    sil_type = silhouette_samples(np.array(embedding_umap), label)
    sil_omic = silhouette_samples(np.array(embedding_umap), tech)
    sil_f1 = (
                2
                * (1 - (sil_omic + 1) / 2)
                * (sil_type + 1)
                / 2
                / (1 - (sil_omic + 1) / 2 + (sil_type + 1) / 2)
            )
    return sil_type.mean(), sil_omic.mean(), sil_f1.mean()

def _average_precision(match: np.ndarray) -> float:
    if np.any(match):
        cummean = np.cumsum(match) / (np.arange(match.size) + 1)
        return cummean[match].mean().item()
    return 0.0

def mean_average_precision(
        x: np.ndarray, y: np.ndarray, neighbor_frac: float = 0.01, **kwargs
) -> float:
    r"""
    Mean average precision

    Parameters
    ----------
    x
        Coordinates
    y
        Cell type labels
    neighbor_frac
        Nearest neighbor fraction
    **kwargs
        Additional keyword arguments are passed to
        :class:`sklearn.neighbors.NearestNeighbors`

    Returns
    -------
    map
        Mean average precision
    """
    k = max(round(y.shape[0] * neighbor_frac), 1)
    nn = NearestNeighbors(
        n_neighbors=min(y.shape[0], k + 1), **kwargs
    ).fit(x)
    nni = nn.kneighbors(x, return_distance=False)
    match = np.equal(y[nni[:, 1:]], np.expand_dims(y, 1))
    return np.apply_along_axis(_average_precision, 1, match).mean().item()

def avg_silhouette_width(x: np.ndarray, y: np.ndarray, **kwargs) -> float:
    r"""
    Cell type average silhouette width

    Parameters
    ----------
    x
        Coordinates
    y
        Cell type labels
    **kwargs
        Additional keyword arguments are passed to
        :func:`sklearn.metrics.silhouette_score`

    Returns
    -------
    asw
        Cell type average silhouette width

    Note
    ----
    Follows the definition in `OpenProblems NeurIPS 2021 competition
    <https://openproblems.bio/neurips_docs/about_tasks/task3_joint_embedding/>`__
    """
    return (silhouette_score(x, y, **kwargs).item() + 1) / 2