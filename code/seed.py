import random
import matplotlib.pyplot as plt
from matplotlib.ticker import StrMethodFormatter
import pandas as pd
import numpy as np
import seaborn as sns
import math

import anndata

if __name__ == "__main__":
    adata = anndata.read_h5ad("../data//Mousecortex_protocols.h5ad")
    genes = adata.var_names
    random.seed(2022)
    sample_genes = random.sample(list(genes), 20)

    ## first get the number of cell types
    celltypes = ['Astrocytes', 'Endothelial', 'Interneuron', 'Microglia', 'Neuron', 'Oligodendrocytes', 'Pericyte',
                 'Polydendrocytes']
    num_celltypes = len(celltypes)

    for gene in sample_genes:

        data = adata[:, gene]
        fig, ax = plt.subplots((num_celltypes) // 3 + 1, 3)  ## draw multi-facet plot
        count = 0

        ## plot all the gene expression values for this gene

        sns.kdeplot(np.log10(data.X.toarray().ravel() + 1), bw=0.5, ax=ax[0, 0]).set_title("All Cells")
        for type in celltypes:
            count += 1
            ## extract gene expression for Astrocytes
            astro_adata = data[data.obs['cell.type'] == type]  ## 1412 cells
            sns.kdeplot(np.log10(data.X.toarray().ravel() + 1), bw=0.5, ax=ax[0, 0]).set_title("All Cells")
            sns.kdeplot(np.log10(astro_adata.X.toarray().ravel() + 1), bw=0.5,
                        ax=ax[(math.floor(count / 3)), (count) % 3]).set_title(
                type)  ##extract out gene expressions for this cell type and plot
        #ax[0].set_title(gene)
        plt.suptitle(gene)
        plt.tight_layout()
        plt.show()
        break
        
