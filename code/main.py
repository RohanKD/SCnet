import anndata
import numpy as np
import tensorflow as tf
from sklearn.preprocessing import OneHotEncoder

## import my written library
from code import _utils

if __name__ == '__main__':
    ## load train_adata and test_adata
    # loads in dataset
    adata = anndata.read_h5ad("../data/Mouse_pFC.h5ad")
    # print("adata")
    #print(adata)
    # print("adata.obs")
    # print(adata.obs)
    # print("adata.n_vars")
    # print(adata.n_vars)
    # print("adata.n_obs")
    # print(adata.n_obs)
    # print("SAMPLE")
    # print(adata.obs['Sample'])
    # set(adata.obs['Sample'])
    # print(adata.obs['cell.type'].value_counts())
    # print(adata.X.shape)
    # print(adata.X)
    # print(np.mean(adata.X) ) ## overall mean
    # print(np.mean(adata.X, axis=0)  )## column mean
    # print(np.mean(adata.X, axis=1)  )## row mean
    # print(np.sum(adata.X) ) ## overall sum
    # print(np.sum(adata.X, axis=0)  )## column sum
    # print(np.sum(adata.X, axis=1) ) ## row sum

    ## you might want to have some visualization
    # plt.hist(adata.X[:, 0])  ## plot the distribution of first gene
    #     # plt.show()
    #     # plt.hist(np.log(adata.X[:, 0] + 0.1))  ## make a log-transformed change and check what has changed
    #     # plt.show()
    # creates train and test
    train_adata = adata[adata.obs['Experiment'] == "Cortex1"]
    #sample2_adata = adata[adata.obs['SampleID'] == "Cortex2"]
    #train_adata = anndata.AnnData.concatenate(*[sample1_adata, sample2_adata],
                                              #join="inner")
    # train_adata = adata[adata.obs['Sample'] == 'PFCSample11']
    test_adata = adata[adata.obs['Experiment'] == 'Cortex2']

    ## extract common genes first
    common_genes = set(train_adata.var_names).intersection(set(test_adata.var_names))
    train_adata = train_adata[:, list(common_genes)]
    test_adata = test_adata[:, list(common_genes)]

    ## preprocess the training data
    train_adata = _utils._process_adata(train_adata, process_type='train')
    train_adata = _utils._select_feature(train_adata,
                                         fs_method='F-test',
                                         num_features=1000)  ## use F-test to select 1000 informative genes
    train_adata = _utils._scale_data(train_adata)  ## center-scale
    _utils._visualize_data(train_adata, output_dir=".",
                           prefix="traindata_vis")  ## visualize cell types with selected features on a low dimension (you might need to change some parameters to let them show all the cell labels)
    ## train an MLP model on it
    MLP_DIMS = _utils.MLP_DIMS  ## get MLP structure from _utils.py
    x_train = _utils._extract_adata(train_adata)
    enc = OneHotEncoder(handle_unknown='ignore')
    y_train = enc.fit_transform(train_adata.obs[[_utils.Celltype_COLUMN]]).toarray()
    mlp = _utils._init_MLP(x_train, y_train, dims=MLP_DIMS,
                           seed=_utils.RANDOM_SEED)
    mlp.compile()
    mlp.fit(x_train, y_train)
    mlp.model.save('./trained_MLP')  ## save the model so that you can load and play with it
    encoders = dict()
    for idx, cat in enumerate(enc.categories_[0]):
        encoders[idx] = cat
    # set(adata.obs['Sample'])
    ## preprocess the test data and predict cell types
    test_adata = _utils._process_adata(test_adata, process_type='test')
    test_adata = test_adata[:,
                 list(train_adata.var_names)]  ## extract out the features selected in the training dataset
    test_data_mat = _utils._extract_adata(test_adata)
    test_data_mat = (test_data_mat - np.array(train_adata.var['mean'])) / np.array(train_adata.var['std'])
    y_pred = tf.nn.softmax(mlp.model.predict(test_data_mat)).numpy()
    pred_celltypes = _utils._prob_to_label(y_pred, encoders)
    test_adata.obs[_utils.PredCelltype_COLUMN] = pred_celltypes
    print(type(pred_celltypes))
    ## let us evaluate the performance --> luckily you will have the accuracy over 99%
    from sklearn.metrics import accuracy_score, adjusted_rand_score, f1_score

    print("Overall Accuracy:",
          accuracy_score(test_adata.obs[_utils.Celltype_COLUMN], test_adata.obs[_utils.PredCelltype_COLUMN]))
    print("ARI:",
          adjusted_rand_score(test_adata.obs[_utils.Celltype_COLUMN], test_adata.obs[_utils.PredCelltype_COLUMN]))
    print("Macro F1:",
          f1_score(test_adata.obs[_utils.Celltype_COLUMN], test_adata.obs[_utils.PredCelltype_COLUMN], average='macro'))
    ## a lot more evaluation metrics can be found on sklearn.metrics and you can explore with them
