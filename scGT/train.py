import torch
import torch.nn as nn
import torch.nn.functional as F
import dataset
from eval import evaluate_cpu, eval_acc
from utils import hard_loss, query_graph_loss, save_model, move_query_to_reference, sub_graph
import matplotlib.pyplot as plt
import random

def model_train(args, dataset, model, split_idx, device, x, n, adjs, adj_loss_inter, adj_loss_intra2):
    l = float('inf') 
    num = 0
    pro = 0
    criterion = nn.NLLLoss()
    eval_func = eval_acc

    rna_idx = split_idx['train']
    atac_idx = torch.arange(dataset.n_data1, n)
    train_atac_idx = torch.empty(0)
    label_train = dataset.label.squeeze(1).clone()
    
    model.reset_parameters()
    optimizer = torch.optim.Adam(model.parameters(),weight_decay=args.weight_decay, lr=args.lr)
    # scheduler = torch.optim.lr_scheduler.MultiStepLR(optimizer, milestones=[200], gamma=0.5)
        
    for epoch in range(args.maxepochs):
        model.to(device)
        model.train()

        if (args.num_batch > 1):
            random_seed = random.randint(0, 10000)
            torch.manual_seed(random_seed)

        num_batch1 = rna_idx.size(0) // args.num_batch
        num_batch2 = atac_idx.size(0) // args.num_batch
        num_batch3 = train_atac_idx.size(0) // args.num_batch

        idx1 = torch.randperm(rna_idx.size(0)) 
        idx2 = torch.randperm(atac_idx.size(0)) 
        idx3 = torch.randperm(train_atac_idx.size(0)) 
    
        label_train = label_train.to(device)    
        L = 0
        L1all = 0
        L2all = 0
        L3all = 0
        
        for i in range(args.num_batch):
            idx_i_rna = rna_idx[idx1[i*num_batch1:(i+1)*num_batch1]]
            idx_i_atac = atac_idx[idx2[i*num_batch2:(i+1)*num_batch2]]
            idx_i_train_atac = train_atac_idx[idx3[i*num_batch3:(i+1)*num_batch3]]
            idx_i = torch.cat((idx_i_rna, idx_i_train_atac, idx_i_atac), dim=0).long()
            x_i = x[idx_i].to(device)

            sub_edge_inter, sub_edge_intra2, adjs_i = sub_graph(args, idx_i, adj_loss_inter, adj_loss_intra2, n, adjs, device)
        
            optimizer.zero_grad()
        
            out_i, _, z1, z2 = model(x_i, adjs_i)
            out_i = F.log_softmax(out_i, dim=1)
            p = F.softmax(out_i, dim=1)
        
            idx_cross_entropy = torch.cat((idx_i_rna, idx_i_train_atac), dim=0).long()
            loss1 = criterion(out_i[0:idx_cross_entropy.shape[0]], label_train[idx_cross_entropy])
            loss2 = hard_loss(sub_edge_inter, z2, args.lamda1)
            loss3 = query_graph_loss(p, sub_edge_intra2, args.lamda2)
        
            loss = loss1 + loss2 + loss3
            L += loss
            L1all += loss1
            L2all += loss2
            L3all += loss3

            loss.backward()
            optimizer.step()
            # scheduler.step()

        if epoch % args.eval_step == 0: #args.eval_step
            result = evaluate_cpu(model, dataset, split_idx, eval_func, criterion, args)

            print(f'Epoch: {epoch:02d}, '
                  f'Loss: {L/args.num_batch:.4f}, '
                  f'Cross Entropy Loss: {L1all/args.num_batch:.4f}, '
                  f'Hard regularity: {L2all/args.num_batch:.4f}, '
                  f'Query Graph regularity: {L3all/args.num_batch:.4f}, '
                  f'Reference: {100 * result:.2f}% ')
            if epoch % 50 == 0 and args.is_move:
                model.eval() 
                model.to(torch.device("cpu"))
                atac_idx, train_atac_idx, label_train, p = move_query_to_reference(model, x, adjs, atac_idx, train_atac_idx, label_train, device)
                
        if l > L/args.num_batch:
            l = L/args.num_batch
            num = 0
            model.eval() 
            model.to(torch.device("cpu"))
            query_acc = save_model(args, dataset, model, x, adjs, epoch) 
            model.to(device)
        else:
            num = num + 1
        
        # pro = sum(torch.max(p[dataset.n_data1:,:], axis=1)[0] > 0.8)/dataset.n_data2
        # print(pro)
        if num >= args.early_stop or torch.isnan(loss):# or pro > 0.8:
            print (f'Accuracy (Query data): {100 * query_acc:.2f}%')
            if torch.isnan(loss):
                print (f'loss is nan!')
            break;
