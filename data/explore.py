### script to explore the data

import anndata
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

## read in 
adata = anndata.read_h5ad("Mouse_pFC.h5ad")
## take a look
adata
adata.n_vars
adata.n_obs

adata.obs  ## look at the metadata for each cell
adata.var  ## look at the metadata for each gene
## Note these are in the format of pandas dataframe 
## You can check this simple tutorial to get familiar with pandas dataframe: https://pandas.pydata.org/docs/user_guide/10min.html
adata.obs['Sample']  ## look at the sample column
set(adata.obs['Sample'])  ## check what samples are in the column

## say we want to subset all cells with Sample == P21Sample1
adata.obs['Sample'] == "P21Sample1"  ## this gives you the indicator about which row to select
adata[adata.obs['Sample'] == "P21Sample1"]  ## this gives you the anndata object which meets the condition

## combine multiple samples with same dimension
sample1_adata = adata[adata.obs['Sample'] == "P21Sample1"]
sample2_adata = adata[adata.obs['Sample'] == "P21Sample2"]
combined_adata = anndata.AnnData.concatenate(*[sample1_adata, sample2_adata],
                                             join="inner")  ## this method would add a column called "batch": https://anndata.readthedocs.io/en/latest/generated/anndata.AnnData.concatenate.html#anndata.AnnData.concatenate
## you can also refer here: https://anndata.readthedocs.io/en/latest/concatenation.html
## for alternative concatenation/merge methods

## combine multiple samples with different dimension
sample1_adata = adata[adata.obs['Sample'] == "P21Sample1", 0:500]  ## with first 500 genes
sample2_adata = adata[adata.obs['Sample'] == "P21Sample2", 0:1000]  ## with first 1000 genes
combined_adata = anndata.AnnData.concatenate(*[sample1_adata, sample2_adata],
                                             join="inner")  ## will automatically find the common genes to concatenate cells

## count cell types
adata.obs['cell.type'].value_counts()  ## count cell types in the whole data
sample1_adata.obs['cell.type'].value_counts()  ## count cell types in one sample
## if you need to count for each sample, you need to do a "for loop"

## look at the content
adata.X  ## this is a numpy array which stores the gene expression matrix, you can refer to: https://numpy.org/doc/stable/user/quickstart.html
adata.X.shape  ## check the shape of the data

## calculate mean, sum or other statistics..
np.mean(adata.X)  ## overall mean
np.mean(adata.X, axis=0)  ## column mean
np.mean(adata.X, axis=1)  ## row mean
np.sum(adata.X)  ## overall sum
np.sum(adata.X, axis=0)  ## column sum
np.sum(adata.X, axis=1)  ## row sum

## you might want to have some visualization
plt.hist(adata.X[:, 0])  ## plot the distribution of first gene
plt.show()
plt.hist(np.log(adata.X[:, 0] + 0.1))  ## make a log-transformed change and check what has changed
plt.show()
