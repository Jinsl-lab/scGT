## scGT: Integration algorithm for single-cell RNA-seq and ATAC-seq based on Graph Transformer

This work proposes scGT, based on Graph Transformer for the integration of scRNA-seq and scATAC-seq data. A robust hybrid graph construction method is included inside.

### Datasets

Human PBMC: [data](https://github.com/SydneyBioX/scJoint/blob/main/data.zip), [Fragment files](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE156478)

Mouse colon: [GSE207308](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE207308)

Mouse brain cortex: [GSE126074](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE126074), SRP18352118

Human myocardial infarction: [data](https://cellxgene.cziscience.com/collections/8191c283-0816-424b-9b61-c3e1d6258a77), peak-by-cell matrix([part1](https://zenodo.org/records/6578553), [part2](https://zenodo.org/record/6578617))

Human fetal atlas: [GSE156793](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE156793), [GSE149683](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE149683)


### How to run our codes?

1. First, run graph_build/main_graph.R

2. Then, you can run scGT in batch_main.ipynb or main.ipynb
