import numpy as np
from sklearn.manifold import TSNE
import matplotlib.pyplot as plt
import torch.nn.functional as F
import torch
import umap

def umap_emb(args, model, dataset, x, adjs):
    model.load_state_dict(torch.load(args.data_dir + 'model/' + f'{args.dataset}.pkl'))
    model.to(torch.device("cpu"))
    model.eval()

    out_i, link_loss_, z1, z2 = model(x, adjs)
    out_i = F.log_softmax(out_i, dim=1)
    p = F.softmax(out_i, dim=1)

    embedding = z2.cpu()
    embedding = embedding[0].detach().numpy()

    label = dataset.label.cpu()
    label = label.squeeze().numpy() 

    tech = dataset.metadata.cpu()
    tech = tech.squeeze().numpy() 

    out = p.cpu() 
    pre = out.argmax(dim=-1, keepdim=True).detach().squeeze().numpy()

    num_celltype = dataset.num_celltype 
    print(num_celltype)

    umap_model = umap.UMAP(n_neighbors=5, n_components=2, metric='euclidean')
    embedding_umap = umap_model.fit_transform(embedding)
    
    return embedding_umap, label, pre, tech