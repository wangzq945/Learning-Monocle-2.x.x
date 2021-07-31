# The input data 

1. 在 `Linux Shell` 中，运行以下代码，下载 `input data`，必要时解压
2. 在浏览器中，输入以下 `url` ，下载 `input data`，并存储在 `data` 文件夹，必要时解压

```sh
mkdir data
```

```sh
# for: monocle
# dataset of Peripheral Blood Mononuclear Cells (PBMC) freely available from 10X Genomics
file="pbmc3k_filtered_gene_bc_matrices.tar.gz"
url="https://cf.10xgenomics.com/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz"
cd data
wget $url
tar -zxvf $file
cd ../
```
