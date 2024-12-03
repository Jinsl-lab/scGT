### PCA版本，其实就是内图只按K值生成一个，不分那么多了 ###

################# 图间连接过滤 #################
# 以图内高准确率的连接结果，借用邻域，对图间连接进行过滤 ，以提高图间连接的准确率

suppressMessages(library(igraph))
suppressMessages(library(dplyr))
suppressMessages(library(Matrix))

filter4_weight <- function(dir) {
  inter <- read.csv(paste0(dir, "/inter_graph0.csv"))
  label1 <- read.csv(paste0(dir, "/label1.csv")) 
  label2 <- read.csv(paste0(dir, "/label2.csv"))
  intra1 <- read.csv(paste0(dir, "/intra_graph1.csv"))
  intra2 <- read.csv(paste0(dir, "/intra_graph2.csv"))
  label1 <- as.matrix(label1)
  label2 <- as.matrix(label2)
  n1 <- length(label1)
  n2 <- length(label2)
  label <- rbind(label1, label2)
  
  # 生成两内图的邻接矩阵
  edges_inter <- as.matrix(inter[, c("V1", "V2")])+1
  edges_intra1 <- as.matrix(intra1[, c("V1", "V2")])
  edges_intra2 <- as.matrix(intra2[, c("V1", "V2")])

  pp <- 1
  
  print(paste("------------------------开始矩阵计算: ", Sys.time())) 
    
  edges_intra1 <- t(as.matrix(intra1[, c("V1", "V2")]))
  edges_intra1 <- list(v1 = edges_intra1[1, ]+1, v2 = edges_intra1[2, ]+1)
  intra1_adj <- graph_from_data_frame(edges_intra1, directed = TRUE)
  intra1_adj <- as_adjacency_matrix(intra1_adj)
  # adj_2hop <- intra1_adj %*% intra1_adj
  # adj2 <- (adj_2hop > 0) | (intra1_adj > 0)
  # adj2 <- as(adj2, "sparseMatrix") - intra1_adj
  # intra1_adj_now <- intra1_adj + adj2
  # adj_2hop <- intra1_adj %*% intra1_adj %*% intra1_adj
  # adj2 <- (adj_2hop > 0) | (intra1_adj_now > 0)
  # adj2 <- as(adj2, "sparseMatrix") - intra1_adj
  # intra1_adj_now <- intra1_adj + adj2
  # intra1_adj <- intra1_adj_now
  
  edges_intra2 <- t(as.matrix(intra2[, c("V1", "V2")]))
  edges_intra2 <- list(v1 = edges_intra2[1, ]+1, v2 = edges_intra2[2, ]+1)
  intra2_adj <- graph_from_data_frame(edges_intra2, directed = TRUE)
  intra2_adj <- as_adjacency_matrix(intra2_adj)
  adj_2hop <- intra2_adj %*% intra2_adj
  adj2 <- (adj_2hop > 0) | (intra2_adj > 0)
  adj2 <- as(adj2, "sparseMatrix") - intra2_adj
  intra2_adj_now <- intra2_adj + adj2*0.9
  adj_2hop <- intra2_adj %*% intra2_adj %*% intra2_adj
  adj2 <- (adj_2hop > 0) | (intra2_adj_now > 0)
  adj2 <- as(adj2, "sparseMatrix") - (intra2_adj_now > 0)
  intra2_adj_now <- intra2_adj_now + adj2*0.81
  adj_2hop <- intra2_adj %*% intra2_adj %*% intra2_adj %*% intra2_adj
  adj2 <- (adj_2hop > 0) | (intra2_adj_now > 0)
  adj2 <- as(adj2, "sparseMatrix") - (intra2_adj_now > 0)
  intra2_adj_now <- intra2_adj_now + adj2*0.729
  intra2_adj <- intra2_adj_now
  # adj_2hop <- intra2_adj %*% intra2_adj %*% intra2_adj %*% intra2_adj %*% intra2_adj
  # adj2 <- (adj_2hop > 0) | (intra2_adj_now > 0)
  # adj2 <- as(adj2, "sparseMatrix") - intra2_adj_now
  # intra2_adj_now <- intra2_adj_now + adj2
  # intra2_adj <- intra2_adj_now
  
  print(paste("------------------------结束矩阵计算: ", Sys.time()))
    
  ### 先考虑data1中每个节点对应的标签在邻域中的占比
  ###（以此衡量它为边界点的可能，边界点所对应的锚点对有更大可能匹配错误）
  
  # 去除邻接矩阵的对角线，即节点自己(这个版本费内存)
  # intra1_adj_undiag <- intra1_adj
  # diag(intra1_adj_undiag) <- 0  
  # # 得到每个节点的邻接节点的类别，并计算节点本身类别在其中的比例。
  # # 若没有邻域，则暂定默认为0.5。有邻域但类型不匹配的为0。
  # p_value <- rep(0.0, n1)
  # indices <- apply(intra1_adj_undiag != 0, 1, which)
  # n_indices <- sapply(indices, length)
  # for (g in 1:n1) {
  #   if (n_indices[g] > 0) {
  #     p_value[g] <- (table(label1[indices[[g]]])[label1[g]]) / n_indices[g]
  #   }
  # }
  # p_value[is.na(p_value)] <- 0
    
    
  diag(intra1_adj) <- 0
  # 得到每个节点的邻接节点的类别，并计算节点本身类别在其中的比例
  # 若没有邻域，则暂定默认为0.5。有邻域但类型不匹配的为0。
  p_value <- rep(0.0, n1)
  indices <- which(intra1_adj != 0, arr.ind = TRUE)
  for (g in 1:n1) {
    adj_nodes <- indices[indices[, 1] == g, 2]
    n_adj <- length(adj_nodes)
    if (n_adj > 0) {
      type_count <- table(label1[adj_nodes])
      p_value[g] <- type_count[label1[g]] / n_adj
    }
  }
  p_value[is.na(p_value)] <- 0  
    
  print(paste("------------------------结束1: ", Sys.time()))        
  
  # 组合得到data1节点的类型与相对于data2锚点来说的可靠比
  data1_p_celltype <- list(celltype =  label1[,1], p_value = p_value)
  
  ### 下面由锚点对得到data2中各节点的预估类别
  ### 将锚点对中data1的类别p_value值加在data2-类别矩阵相应的位置上
  ### 最后取最大值对应的类别作为data2的类别，若都为0，则类别为“UN”
  
  # 定义data2-类别矩阵
  data2_celltype_matrix <- matrix(0, nrow = n2, ncol = length(table(label1)), 
                                  dimnames = list(c(1:n2), names(table(label1))))
  for (k in 1:nrow(edges_inter)) {
    j1 <- edges_inter[k, 1]
    j2 <- edges_inter[k, 2] - n1
    data2_celltype_matrix[j2, data1_p_celltype$celltype[j1]] <- data2_celltype_matrix[j2, data1_p_celltype$celltype[j1]] + data1_p_celltype$p_value[j1]
  }

  print(paste("------------------------结束2: ", Sys.time())) 
  
  # 下面确定data2的预估类别
  data2_celltype_pre <- c()
  for (i in 1:n2) {
    positive <- which(data2_celltype_matrix[i,] > 0)
    index <- which.max(positive)
    if (length(index) != 0) {
      data2_celltype_pre <- c(data2_celltype_pre, names(index)[1])
    } else {
      data2_celltype_pre <- c(data2_celltype_pre, NA)
    }
  }
  
  print("预估的label2与实际值的对比：")
  print(table(data2_celltype_pre == label2))

  print(paste("------------------------data2预估类别结束: ", Sys.time()))    
    
  ### 下面开始过滤流程
  ####
  max_cell_types <- sapply(1:nrow(intra2_adj), function(node2) {
    index <- which(intra2_adj[node2, ] != 0)
    b <- data2_celltype_pre[index]
    weight <- intra2_adj[node2, ][index]
    if (node2 %% 1000 == 0) {
        print(node2)
    }
    if (length(b) > 0) {
      sum_weight <- tapply(weight, b, sum)
      max_indices <- which(sum_weight == max(sum_weight))
      max_cell_type <- names(sum_weight)[max_indices]
      max_cell_type
    } else {
      NULL
    }
  })

  key <- logical(nrow(edges_inter))
  for (i in seq_len(nrow(edges_inter))) {
    node1 <- edges_inter[i, 1]
    node2 <- edges_inter[i, 2] - n1
    ca <- label1[node1]
    cb <- max_cell_types[[node2]]  
    key[i] <- (!is.null(cb) && (ca %in% cb))
  }
  # 计算相等
  a <- label[inter$V1 + 1]
  b <- label[inter$V2 + 1]
  before <- a == b    
  
  print(paste("------------------------结束3: ", Sys.time())) 
  
  after <- c()
  index <- c()
  for (i in 1:length(before)) {
    if (key[i] == TRUE) {
      after <- c(after, before[i])
      index <- c(index, i)
    }
  }
  print("before table: ")
  print(table(before))
  print("after table: ")
  print(table(after))
  
  print("before acc: ")
  print(sum(before == TRUE)/length(before))
  print("after acc: ")
  print(sum(after == TRUE)/length(after))

  print(paste("------------------------过滤结束: ", Sys.time()))
    
  # 更新图间链接
  edges_inter_after <- edges_inter[index,]-1
  edges_inter_after <- list(V1 = edges_inter_after[,1], V2 = edges_inter_after[,2])
  write.csv(edges_inter_after, paste0(dir, "/inter_graph.csv"))
}

filter4_notweight <- function(dir) {
  inter <- read.csv(paste0(dir, "/inter_graph0.csv"))
  label1 <- read.csv(paste0(dir, "/label1.csv")) 
  label2 <- read.csv(paste0(dir, "/label2.csv"))
  intra1 <- read.csv(paste0(dir, "/intra_graph1.csv"))
  intra2 <- read.csv(paste0(dir, "/intra_graph2.csv"))
  label1 <- as.matrix(label1)
  label2 <- as.matrix(label2)
  n1 <- length(label1)
  n2 <- length(label2)
  label <- rbind(label1, label2)
  
  # 生成两内图的邻接矩阵
  edges_inter <- as.matrix(inter[, c("V1", "V2")])+1
  edges_intra1 <- as.matrix(intra1[, c("V1", "V2")])
  edges_intra2 <- as.matrix(intra2[, c("V1", "V2")])

  pp <- 1
  
  print(paste("------------------------开始矩阵计算: ", Sys.time())) 
    
  edges_intra1 <- t(as.matrix(intra1[, c("V1", "V2")]))
  edges_intra1 <- list(v1 = edges_intra1[1, ]+1, v2 = edges_intra1[2, ]+1)
  intra1_adj <- graph_from_data_frame(edges_intra1, directed = TRUE)
  intra1_adj <- as_adjacency_matrix(intra1_adj)
  
  edges_intra2 <- t(as.matrix(intra2[, c("V1", "V2")]))
  edges_intra2 <- list(v1 = edges_intra2[1, ]+1, v2 = edges_intra2[2, ]+1)
  intra2_adj <- graph_from_data_frame(edges_intra2, directed = TRUE)
  intra2_adj <- as_adjacency_matrix(intra2_adj)
  adj_2hop <- intra2_adj %*% intra2_adj %*% intra2_adj %*% intra2_adj# %*% intra2_adj %*% intra2_adj
  adj2 <- (adj_2hop > 0) * 1
  intra2_adj <- as(adj2, "sparseMatrix")
  
  print(paste("------------------------结束矩阵计算: ", Sys.time()))
    
  ### 先考虑data1中每个节点对应的标签在邻域中的占比
  ###（以此衡量它为边界点的可能，边界点所对应的锚点对有更大可能匹配错误）
  
  # 去除邻接矩阵的对角线，即节点自己
  intra1_adj_undiag <- intra1_adj
  diag(intra1_adj_undiag) <- 0
  
  # 得到每个节点的邻接节点的类别，并计算节点本身类别在其中的比例。
  # 若没有邻域，则暂定默认为0.5。有邻域但类型不匹配的为0。
  p_value <- rep(0.0, n1)
  indices <- apply(intra1_adj_undiag != 0, 1, which)
  n_indices <- sapply(indices, length)
  for (g in 1:n1) {
    if (n_indices[g] > 0) {
      p_value[g] <- (table(label1[indices[[g]]])[label1[g]]) / n_indices[g]
    }
  }
  p_value[is.na(p_value)] <- 0

  print(paste("------------------------结束1: ", Sys.time()))        
  
  # 组合得到data1节点的类型与相对于data2锚点来说的可靠比
  data1_p_celltype <- list(celltype =  label1[,1], p_value = p_value)
  
  ### 下面由锚点对得到data2中各节点的预估类别
  ### 将锚点对中data1的类别p_value值加在data2-类别矩阵相应的位置上
  ### 最后取最大值对应的类别作为data2的类别，若都为0，则类别为“UN”
  
  # 定义data2-类别矩阵
  data2_celltype_matrix <- matrix(0, nrow = n2, ncol = length(table(label1)), 
                                  dimnames = list(c(1:n2), names(table(label1))))
  for (k in 1:nrow(edges_inter)) {
    j1 <- edges_inter[k, 1]
    j2 <- edges_inter[k, 2] - n1
    data2_celltype_matrix[j2, data1_p_celltype$celltype[j1]] <- data2_celltype_matrix[j2, data1_p_celltype$celltype[j1]] + data1_p_celltype$p_value[j1]
  }

  print(paste("------------------------结束2: ", Sys.time())) 
  
  # 下面确定data2的预估类别
  data2_celltype_pre <- c()
  for (i in 1:n2) {
    positive <- which(data2_celltype_matrix[i,] > 0)
    index <- which.max(positive)
    if (length(index) != 0) {
      data2_celltype_pre <- c(data2_celltype_pre, names(index)[1])
    } else {
      data2_celltype_pre <- c(data2_celltype_pre, NA)
    }
  }
  
  print("预估的label2与实际值的对比：")
  print(table(data2_celltype_pre == label2))

  print(paste("------------------------data2预估类别结束: ", Sys.time()))    
    
  ### 下面开始过滤流程
  ## 不加权的方法
  max_cell_types <- sapply(1:nrow(intra2_adj), function(node2) {
    b <- data2_celltype_pre[which(intra2_adj[node2, ] != 0)]
    if (node2 %% 1000 == 0) {
      print(node2)
    } 
    if (length(b) > 0) {
      table_b <- table(b)
      max_value <- max(table_b)
      names(table_b[table_b == max_value])
    } else {
      NULL
    }
  })
  key <- logical(nrow(edges_inter))
  for (i in seq_len(nrow(edges_inter))) {
    node1 <- edges_inter[i, 1]
    node2 <- edges_inter[i, 2] - n1
    ca <- label1[node1]
    cb <- max_cell_types[[node2]]  
    key[i] <- (!is.null(cb) && (ca %in% cb))
  }
  # 计算相等
  a <- label[inter$V1 + 1]
  b <- label[inter$V2 + 1]
  before <- a == b    
  
  print(paste("------------------------结束3: ", Sys.time())) 
  
  after <- c()
  index <- c()
  for (i in 1:length(before)) {
    if (key[i] == TRUE) {
      after <- c(after, before[i])
      index <- c(index, i)
    }
  }
  print("before table: ")
  print(table(before))
  print("after table: ")
  print(table(after))
  
  print("before acc: ")
  print(sum(before == TRUE)/length(before))
  print("after acc: ")
  print(sum(after == TRUE)/length(after))

  print(paste("------------------------过滤结束: ", Sys.time()))
    
  # 更新图间链接
  edges_inter_after <- edges_inter[index,]-1
  edges_inter_after <- list(V1 = edges_inter_after[,1], V2 = edges_inter_after[,2])
  write.csv(edges_inter_after, paste0(dir, "/inter_graph.csv"))
}


fix_index <- function(dir) {
  inter <- read.csv(paste0(dir, "/inter_graph0.csv"))
  label1 <- read.csv(paste0(dir, "/label1.csv")) 
  intra2 <- read.csv(paste0(dir, "/intra_graph2.csv"))
  
  inter <- as.matrix(inter[, c("V1", "V2")])
  intra2 <- as.matrix(intra2[, c("V1", "V2")])
  
  label1 <- as.matrix(label1)
  n1 <- length(label1)
  
  inter[, 2] <- inter[, 2] + n1
  intra2 <- intra2 + n1
  
  inter <- list(V1 = inter[,1], V2 = inter[,2])
  intra2 <- list(V1 = intra2[,1], V2 = intra2[,2])
  
  write.csv(inter, paste0(dir, "/inter_graph0.csv"))
  write.csv(intra2, paste0(dir, "/intra_graph2.csv"))
}

acc <- function(dir) {
  inter <- read.csv(paste0(dir, "/inter_graph0.csv"))
  label1 <- read.csv(paste0(dir, "/label1.csv")) 
  label2 <- read.csv(paste0(dir, "/label2.csv"))
  intra1 <- read.csv(paste0(dir, "/intra_graph1.csv"))
  intra2 <- read.csv(paste0(dir, "/intra_graph2.csv"))
  label1 <- as.matrix(label1)
  label2 <- as.matrix(label2)
  label <- rbind(label1, label2)
  
  print(paste("inter_graph0 acc: ", caculate(label[inter$V1+1], label[inter$V2+1], 0), "."))
  print(paste("intra_graph1 acc: ", caculate(label[intra1$V1+1], label[intra1$V2+1], length(label1)), "."))
  print(paste("intra_graph2 acc: ", caculate(label[intra2$V1+1], label[intra2$V2+1], length(label2)), "."))
}

# 已经产生inter_graph的用这个
acc2 <- function(dir) {
  inter0 <- read.csv(paste0(dir, "/inter_graph0.csv"))
  inter <- read.csv(paste0(dir, "/inter_graph.csv"))
  label1 <- read.csv(paste0(dir, "/label1.csv")) 
  label2 <- read.csv(paste0(dir, "/label2.csv"))
  intra1 <- read.csv(paste0(dir, "/intra_graph1.csv"))
  intra2 <- read.csv(paste0(dir, "/intra_graph2.csv"))
  label1 <- as.matrix(label1)
  label2 <- as.matrix(label2)
  label <- rbind(label1, label2)
  
  print(paste("inter_graph acc: ", caculate(label[inter$V1+1], label[inter$V2+1], 0), "."))
  print(paste("inter_graph0 acc: ", caculate(label[inter0$V1+1], label[inter0$V2+1], 0), "."))
  print(paste("intra_graph1 acc: ", caculate(label[intra1$V1+1], label[intra1$V2+1], length(label1)), "."))
  print(paste("intra_graph2 acc: ", caculate(label[intra2$V1+1], label[intra2$V2+1], length(label2)), "."))
}

caculate <- function(a, b, n) {
  # a, b分别是连边的两端节点对应的类别
  return((sum(a == b)-n)/(length(a)-n))
}

max_ <- function(celltype, index, p) {
  k <- rep(0,length(celltype))
  names(k) <- celltype
  for (i in length(index)) {
    if (!is.na(index[i]) & !is.na(p[i]))
      k[index[i]] = k[index[i]] + p[i]
  }
  if (sum(k) == 0) return (NULL)
  else return (names(which.max(k)))
}

adj_nhop <- function(adj, weights, n = 4) {
  # 默认为4跳
  adj_hop <- adj * weights[1]
  b <- adj
  for (i in 2:n) {
    b <- b %*% adj
    adj2 <- (b > 0) * 1
    adj2 <- as(adj2, "sparseMatrix")
    adj_hop <- adj_hop + adj2 * weights[i]
  }
  return(adj_hop)
}
