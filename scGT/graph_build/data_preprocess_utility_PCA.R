suppressMessages(library(Seurat)); suppressMessages(library(entropy)); suppressMessages(library(Matrix))
library(anndata)
library(SeuratDisk)


findMNN <- function(knn) {
  n <- nrow(knn)
  d <- ncol(knn)
  graph <- matrix(, nrow = n * d, ncol = 2)
  idx <- 1
  for (i in 1:n) {
    for (j in 1:d) {
      if (i %in% knn[knn[i, j], ]) {
        graph[idx, ] <- c(i, knn[i, j])
        idx <- idx + 1
      }
    }
  }
  graph <- graph[1:(idx-1), ]
  return(graph)
}

GenerateGraph <- function(Dat1,Dat2,Lab1,dir,K_inter,K_intra,n_gene=n_gene,dim=dim){
  object1 <- CreateSeuratObject(counts=Dat1,project = "1",assay = "Data1",
                                min.cells = 0,min.features = 0,
                                names.field = 1,names.delim = "_")
  
  object2 <- CreateSeuratObject(counts=Dat2,project = "2",assay = "Data2",
                                min.cells = 0,min.features =0,names.field = 1,
                                names.delim = "_")
  print(paste("------------------------ Start preprocessing: ", Sys.time()))
  objects <- list(object1,object2)    
  objects1 <- lapply(objects,function(obj){
    obj <- NormalizeData(obj,verbose=F)
    obj <- FindVariableFeatures(obj,
                                selection.method = "vst",
                                nfeatures = n_gene)
    print(paste("------------------------ finish one: ", Sys.time()))
    return(obj)})
    
  features <- SelectIntegrationFeatures(object.list = objects1, nfeatures = n_gene)
  objects1 <- lapply(X = objects1, FUN = function(x) {
    x <- ScaleData(x, features = features)
    x <- RunPCA(x, features = features)
  })

  print(paste("------------------------ Preprocessing completed: ", Sys.time()))

  print(paste("------------------------ Start calculating inter-dataset graph connections: ", Sys.time()))
  #'  Inter-data graph
  object.nn <- FindIntegrationAnchors(object.list = objects1,k.anchor=K_inter, anchor.features = features, reduction = "rpca")
  arc=object.nn@anchors
  d1.arc1=cbind(arc[arc[,4]==1,1],arc[arc[,4]==1,2],arc[arc[,4]==1,3])
  grp1=d1.arc1[d1.arc1[,3]>0,1:2]-1

  objects1 <- lapply(X = objects1, FUN = function(obj) {
    obj <- ScaleData(obj,features=rownames(obj))
    obj <- RunPCA(obj,features=rownames(obj))
  })
  
  #'  Intra-data graph  
  print(paste("------------------------ Start calculating intra-dataset graph connections: ", Sys.time()))
  d1 <- FindNeighbors(objects1[[1]], dims = 2:dim, k.param = K_intra, return.neighbor = TRUE) 
  d1.knn <- d1@neighbors$Data1.nn@nn.idx
  d1.grp <- findMNN(d1.knn)-1
  
  d2 <- FindNeighbors(objects1[[2]], dims = 2:dim, k.param = K_intra, return.neighbor = TRUE) 
  d2.knn <- d2@neighbors$Data2.nn@nn.idx
  d2.grp <- findMNN(d2.knn)-1

  final <- list(inteG=grp1,intraG1=d1.grp,intraG2=d2.grp)
  return (final)
}

save_processed_data <- function(count.list,label.list,dir="input",K_inter=30,K_intra=10,n_gene=n_gene,dim=50){
  print(paste("------------------------begin: ", Sys.time()))
  dir.create(dir);
  rna <- CreateSeuratObject(count.list[[1]], assay = "RNA")
  atac <- CreateSeuratObject(count.list[[2]], assay = "atac")
  rna$cell_type <- label.list[[1]]
  atac$cell_type <- label.list[[2]]

  rna[["RNA"]] <- as(object = rna[["RNA"]], Class = "Assay")
  SaveH5Seurat(rna, filename = paste0(dir,'/rna.h5Seurat'))
  Convert(paste0(dir,'/rna.h5Seurat'), dest = "h5ad")

  atac[["atac"]] <- as(object = atac[["atac"]], Class = "Assay")
  SaveH5Seurat(atac, filename = paste0(dir,'/atac.h5Seurat'))
  Convert(paste0(dir,'/atac.h5Seurat'), dest = "h5ad")
  rm(rna, atac)
  gc()
  print(paste("------------------------ Data conversion ended: ", Sys.time()))
  
  new.dat1 <- count.list[[1]]; new.dat2 <- count.list[[2]]
  new.lab1 <- label.list[[1]]; new.lab2 <- label.list[[2]]
  graphs <- suppressWarnings(GenerateGraph(Dat1=new.dat1,Dat2=new.dat2,
                                            Lab1=new.lab1,K_inter=K_inter,K_intra=K_intra,n_gene=n_gene,dim=dim,
                                            dir = dir))
  write.csv(graphs[[1]],file=paste0(dir,'/inter_graph0.csv'),quote=F,row.names=T)
  write.csv(graphs[[2]],file=paste0(dir,'/intra_graph1.csv'),quote=F,row.names=T)
  write.csv(graphs[[3]],file=paste0(dir,'/intra_graph2.csv'),quote=F,row.names=T)
  write.csv(new.lab1,file=paste0(dir,'/label1.csv'),quote=F,row.names=F)
  write.csv(new.lab2,file=paste0(dir,'/label2.csv'),quote=F,row.names=F)
}
