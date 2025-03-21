suppressMessages(library(igraph))
suppressMessages(library(dplyr))
suppressMessages(library(Matrix))
library(foreach)
library(doParallel)


filter4_weight_linux <- function(dir, n_ = 12) {
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
  
  edges_inter <- as.matrix(inter[, c("V1", "V2")])+1
  edges_intra1 <- as.matrix(intra1[, c("V1", "V2")])
  edges_intra2 <- as.matrix(intra2[, c("V1", "V2")])
    
  edges_intra1 <- t(as.matrix(intra1[, c("V1", "V2")]))
  edges_intra1 <- list(v1 = edges_intra1[1, ]+1, v2 = edges_intra1[2, ]+1)
  intra1_adj <- graph_from_data_frame(edges_intra1, directed = TRUE)
  intra1_adj <- as_adjacency_matrix(intra1_adj)
  
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
    
  cl <- makeCluster(n_) 
  registerDoParallel(cl)
  p_value <- rep(0.0, n1)
  indices <- which(intra1_adj != 0, arr.ind = TRUE)
    
  # clusterExport(cl, c('indices', 'label1'))
  Sys.sleep(5)    

  log_file <- "progress_log.txt"
  cat("Progress1 started at: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n", file = log_file, append = TRUE)
  p_value <- foreach(g = 1:n1, .combine = c, .packages = 'base', .export = c('indices', 'label1')) %dopar% {
    adj_nodes <- indices[indices[, 1] == g, 2]  
    n_adj <- length(adj_nodes)  
    p <- 0.0  
    if (n_adj > 0) {
      type_count <- table(label1[adj_nodes])  
      p <- type_count[label1[g]] / n_adj  
    }
    if (g %% 1000 == 0) {
      cat("Processed node: ", g, " at ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n", file = log_file, append = TRUE) 
      # flush(stdout()) 
    }
    p  
  }
  p_value[is.na(p_value)] <- 0
  stopCluster(cl)  
  
  saveRDS(p_value, "p_value.rds")
    
  data1_p_celltype <- list(celltype =  label1[,1], p_value = p_value)
  
  data2_celltype_matrix <- matrix(0, nrow = n2, ncol = length(table(label1)), 
                                  dimnames = list(c(1:n2), names(table(label1))))
  for (k in 1:nrow(edges_inter)) {
    j1 <- edges_inter[k, 1]
    j2 <- edges_inter[k, 2] - n1
    data2_celltype_matrix[j2, data1_p_celltype$celltype[j1]] <- data2_celltype_matrix[j2, data1_p_celltype$celltype[j1]] + data1_p_celltype$p_value[j1]
  }

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
  

  cat("Progress2 started at: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n", file = log_file, append = TRUE)

  max_cell_types <- mclapply(1:nrow(intra2_adj), function(node2) {
    tryCatch({
      index <- which(intra2_adj[node2, ] != 0)
      b <- data2_celltype_pre[index]
      weight <- intra2_adj[node2, ][index]
      if (node2 %% 1000 == 0) {
        cat("Processed node: ", node2, " at ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n", file = log_file, append = TRUE)
      }
      if (length(b) > 0) {
        sum_weight <- tapply(weight, b, sum)
        max_indices <- which(sum_weight == max(sum_weight))
        max_cell_type <- names(sum_weight)[max_indices]
        max_cell_type
      } else {
        NULL
      }
    }, error = function(e) {
      cat("Error in task for node ", node2, ": ", e$message, "\n", file = log_file, append = TRUE)
      NULL 
    })
  }, mc.cores = n_)  

  saveRDS(max_cell_types, "max_cell_types.rds")
    
  key <- logical(nrow(edges_inter))
  for (i in seq_len(nrow(edges_inter))) {
    node1 <- edges_inter[i, 1]
    node2 <- edges_inter[i, 2] - n1
    ca <- label1[node1]
    cb <- max_cell_types[[node2]]  
    key[i] <- (!is.null(cb) && (ca %in% cb))
  }
  
  index <- c()
  for (i in 1:length(inter$V1)) {
    if (key[i] == TRUE) {
      index <- c(index, i)
    }
  }
    
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

acc_intra_graph1 <- function(dir) {
  label1 <- read.csv(paste0(dir, "/label1.csv")) 
  intra1 <- read.csv(paste0(dir, "/intra_graph1.csv"))
  label1 <- as.matrix(label1)
  
  print(paste("intra_graph1 acc: ", caculate(label1[intra1$V1+1], label1[intra1$V2+1], length(label1)), "."))
}

caculate <- function(a, b, n) {
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