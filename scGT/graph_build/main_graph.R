# nohup Rscript 图链接生成.R>32w_graph_make.log 2>&1 &   # 后台运行

# note that
# count.list list of reference (scRNA-seq) data and query (scATAC-seq) data; rows are genes and columns are cells. And they have the same genes.
# label.list list of reference label and query label, both are data frames with rownames identical with colnames of data; the first column is cell type



library(future)
source("data_preprocess_utility_PCA.R")
source("filter_PCA.R")
# source("filter_PCA_linux.R")

plan(strategy = "multicore")
options(future.globals.maxSize = 100000 * 1024^14)

count <- readRDS("count.list.rds")
label <- readRDS("label.list.rds")

# save_processed_data(count, label, Rgraph=TRUE, dir = "input1113")

print(paste("------------------------开始修正idx: ", Sys.time()))
fix_index("input")

print(paste("------------------------过滤前准确率: ", Sys.time()))
acc("input")

print(paste("------------------------开始过滤图连接: ", Sys.time()))
filter4_weight("input")

print(paste("------------------------过滤后准确率: ", Sys.time()))
acc2("input")
