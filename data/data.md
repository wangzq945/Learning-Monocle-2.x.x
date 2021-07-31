# The input data 

1. 在 `Linux Shell` 中，运行以下代码，下载 `input data`，必要时解压
2. 在浏览器中，输入以下 `url` ，下载 `input data`，并存储在 `data` 文件夹，必要时解压

```sh
mkdir data
```

* C. elegans dataset (from [Packer & Zhu et al](http://dx.doi.org/10.1101/565549)): 
  * Their study includes a time series analysis of whole developing embyros. 
  * We will examine a small subset of the data which includes most of the neurons.

```sh
# for: monocle
# C. elegans dataset
file1="packer_embryo_expression.rds"
url1="http://staff.washington.edu/hpliner/data/packer_embryo_expression.rds"
file2="packer_embryo_colData.rds"
url2="http://staff.washington.edu/hpliner/data/packer_embryo_colData.rds"
file3="packer_embryo_rowData.rds"
url3="http://staff.washington.edu/hpliner/data/packer_embryo_rowData.rds"
mkdir -p data/packer_embryo
cd data/packer_embryo
wget $url1
wget $url2
wget $url3
cd ../../
```
