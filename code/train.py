import anndata
import os
import numpy as np

import tensorflow as tf
from sklearn.preprocessing import OneHotEncoder
import math
## import my written library
from utils import _utils

if __name__ == '__main__':
    ## load train_adata and test_adata
    # loads in dataset
    adata = anndata.read_h5ad("../data/Mouse_pFC.h5ad") # test
    bdata = anndata.read_h5ad("../data/Mouse_wholebrain_FC.h5ad") # train
    cdata = anndata.read_h5ad("../data/Mousecortex_protocols.h5ad") # train


    samples = []
    for i in adata.obs['Sample']:
        if i not in samples:
            samples.append(i)


    train_num = len(samples) - (math.ceil(len(samples)/10))
    counter = 0
    for sample in samples:
        if counter == 0:
            train_adata = adata[adata.obs['Sample'] == sample]
        elif(counter < train_num):
            train_adata =  anndata.AnnData.concatenate(*[train_adata,adata[adata.obs['Sample'] == sample] ],join="inner")
        elif(counter == train_num):
            test_adata = adata[adata.obs['Sample'] == sample]
        else:
            test_adata = anndata.AnnData.concatenate(*[test_adata,adata[adata.obs['Sample'] == sample] ],join="inner")

        counter += 1
    b_samples = []
    for i in bdata.obs['sampleID']:
        if i not in b_samples:
            b_samples.append(i)

    b_train_num = len(b_samples) - (math.ceil(len(b_samples) / 10))//2
    counter = 0

    b_train_adata = bdata[bdata.obs['sampleID'] == b_samples[0]]
    for sample in b_samples:
        if counter == 0:
            b_train_adata = bdata[bdata.obs['sampleID'] == sample]
        elif (counter < b_train_num):
            b_train_adata = anndata.AnnData.concatenate(*[b_train_adata, bdata[bdata.obs['sampleID'] == sample]],
                                                      join="inner")
        elif (counter == b_train_num):
            b_test_adata = bdata[bdata.obs['sampleID'] == sample]
        elif (counter < b_train_num):
            b_test_adata = anndata.AnnData.concatenate(*[b_test_adata, bdata[bdata.obs['sampleID'] == sample]], join="inner")

        counter += 1
    c_train_adata = cdata[cdata.obs['Experiment']=="Cortex1"]
    c_test_adata = cdata[cdata.obs['Experiment'] == "Cortex2"]
    c_adata = anndata.AnnData.concatenate(*[c_test_adata, c_test_adata], join="inner")

    combined_train_adata = anndata.AnnData.concatenate(*[b_train_adata,b_test_adata], join="inner")
    combined_train_adata = anndata.AnnData.concatenate(*[c_adata, combined_train_adata], join="inner")
    combined_test_adata = anndata.AnnData.concatenate(*[train_adata, test_adata], join="inner")
    common_genes = set(combined_train_adata.var_names).intersection(set(combined_test_adata.var_names))
    combined_train_adata = combined_train_adata[:, list(common_genes)]
    combined_test_adata = combined_test_adata[:, list(common_genes)]
    combined_train_adata = _utils._process_adata(combined_train_adata, process_type='train')
    combined_train_adata = _utils._select_feature(combined_train_adata,
                                         fs_method='F-test',
                                         num_features=1000)  ## use F-test to select 1000 informative genes
    combined_train_adata = _utils._scale_data(combined_train_adata)  ## center-scale
    # _utils._visualize_data(combined_train_adata, output_dir=".",
    #                        prefix="traindata_vis")  ## visualize cell types with selected features on a low dimension (you might need to change some parameters to let them show all the cell labels)
    # train an MLP model on it
    MLP_DIMS = _utils.MLP_DIMS  ## get MLP structure from _utils.py
    x_train = _utils._extract_adata(combined_train_adata)
    enc = OneHotEncoder(handle_unknown='ignore')
    y_train = enc.fit_transform(combined_train_adata.obs[[_utils.Celltype_COLUMN]]).toarray()
    mlp = _utils._init_MLP(x_train, y_train, dims=MLP_DIMS,
                           seed=_utils.RANDOM_SEED)
    mlp.compile()
    mlp.fit(x_train, y_train)
    mlp.model.save('./trained_MLP')  ## save the model so that you can load and play with it
    encoders = dict()
    for idx, cat in enumerate(enc.categories_[0]):
        encoders[idx] = cat
    print("encoders")
    print(encoders)
    # set(adata.obs['Sample'])
    ## preprocess the test data and predict cell types
    combined_test_adata = _utils._process_adata(combined_test_adata, process_type='test')
    combined_test_adata = combined_test_adata[:,
                 list(combined_train_adata.var_names)]  ## extract out the features selected in the training dataset
    test_data_mat = _utils._extract_adata(combined_test_adata)
    test_data_mat = (test_data_mat - np.array(combined_train_adata.var['mean'])) / np.array(combined_train_adata.var['std'])
    y_pred = tf.nn.softmax(mlp.model.predict(test_data_mat)).numpy()
    pred_celltypes = _utils._prob_to_label(y_pred, encoders)
    combined_test_adata.obs[_utils.PredCelltype_COLUMN] = pred_celltypes

    ## let us evaluate the performance --> luckily you will have the accuracy over 99%
    from sklearn.metrics import accuracy_score, adjusted_rand_score, f1_score

    print("Overall Accuracy:",
          accuracy_score(combined_test_adata.obs[_utils.Celltype_COLUMN], combined_test_adata.obs[
              _utils.PredCelltype_COLUMN]))
    print("ARI:",
          adjusted_rand_score(combined_test_adata.obs[_utils.Celltype_COLUMN], combined_test_adata.obs[
              _utils.PredCelltype_COLUMN]))
    print("Macro F1:",
          f1_score(combined_test_adata.obs[_utils.Celltype_COLUMN], combined_test_adata.obs[_utils.PredCelltype_COLUMN], average='macro'))
    ## a lot more evaluation metrics can be found on sklearn.metrics and you can explore with them
    #save feature file
    model_save_dir = "C:\\Users\\rohan\\PycharmProjects\\SciFaitEmory\\code\\models"
    combined_train_adata.var.loc[:, ['mean', 'std']].to_csv("C:\\Users\\rohan\\PycharmProjects\\SciFaitEmory\\code\\models\\features.txt")
    print(combined_train_adata.var)
    ## save enc information
    with open(model_save_dir + os.sep + "onehot_encoder.txt", 'w') as f:
        for idx, cat in enumerate(enc.categories_[0]):
            f.write('%d:%s\n' % (idx, cat))