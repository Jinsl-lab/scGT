library(future)
source("data_preprocess_utility_PCA.R")
source("filter_PCA.R")
# source("filter_PCA_linux.R")  # only linux can use

plan(strategy = "multicore")
options(future.globals.maxSize = 100000 * 1024^14)

count <- readRDS("../../data/human40k/count.list.rds")
label <- readRDS("../../data/human40k/label.list.rds")

save_processed_data(count, label, dir = "../../data/human40k/input", K_inter = 30, K_intra = 10,n_gene=5000, dim = 50)
### 
# dir: The address where the output file is saved.
# K_inter: The number of k in MNN when the inter-datasets graph connections are generated. It is best to get a number of connections similar to the number of cells in the reference or query data.
# K_intra: The number of k in MNN when the intra-datasets graph connections are generated. Usually it's 10.
# n_gene: The number of hypervariable genes selected when constructing the inter-datasets graph connections.
# dim: The number of principal components when constructing the intra-datasets graph connections.
###

fix_index("input")

acc_intra_graph1("input") # Accuracy of the intra-datasets graph connections in reference data.

filter4_weight("input") # filter

