'''
Functions related to predict cell types
'''
import os, sys

import tensorflow as tf
import numpy as np
import pandas as pd
import anndata
## import my package
import _utils_new


## get the logger
# import logging
# logger = logging.getLogger(__name__)


# change args
def predict():
    # find how to load properly
    model = tf.keras.models.load_model("my_model")

    # load feature file and onehotencoder enumeration file, find them first by training
    # CHANGE
    feature_file = "C:\\Users\\rohan\\PycharmProjects\\SciFaitEmory\\code\\models\\features.txt"
    encoder_file = "C:\\Users\\rohan\\PycharmProjects\\SciFaitEmory\\code\\models\\encoders.txt"

    features = pd.read_csv(feature_file, sep='\t', header=0, index_col=0)
    encoders = {}
    with open(encoder_file) as f:
        for line in f:
            line_info = line.strip().split(':')
            encoders[int(line_info[0])] = line_info[1]

    ## load input data
    print("Loading data... \n This may take a while depending on your data size..")
    # load in from webapp, add parameter for filepath
    target_data = anndata.read_h5ad("../data/Mouse_pFC.h5ad")
    # change this file from args
    # if '.csv' in args.input:
    #     test_adata = _utils._csv_data_loader(args.input)
    # else:
    test_adata = _utils_new._COOmtx_data_loader("../data/Mouse_pFC.h5ad")
    ## process test adata
    test_adata = _utils_new._process_adata(test_adata, process_type='test')

    ## fill in the data with the same order of features
    feature_idx = []
    NA_idx = []
    for f_idx, feature in enumerate(features.index):
        find_flag = False
        for test_idx, gene in enumerate(test_adata.var_names):
            if gene == feature:
                feature_idx.append(test_idx)
                find_flag = True
                break
        if not find_flag:
            feature_idx.append(-1)
            NA_idx.append(f_idx)
    print("%d genes from reference data are found in target.\n" % (len(features) - len(NA_idx)))

    if len(NA_idx) > 0.1 * len(features):
        print("Warnings: too few genes found in target and this will result in inaccurate prediction.")
    if -1 in feature_idx:
        print("Warnings: since some feature does not exist in target dataset. We will fill in 0s for those columns.")
        ## first replace those unique genes with index
        curated_feature_idx = np.array(feature_idx)
        curated_feature_idx[NA_idx] = 0
        test_adata = test_adata[:, curated_feature_idx].copy()
        test_adata.var_names.values[NA_idx] = ["GenesNotFound-" + str(i) for i, NA_item in
                                               enumerate(NA_idx)]  ## change gene names
        test_adata_X = test_adata.X
        test_adata_X[:, NA_idx] = 0
        test_adata.X = test_adata_X
    else:
        test_adata = test_adata[:, feature_idx]
    print("Data shape after processing: %d cells X %d genes" % (test_adata.shape[0], test_adata.shape[1]))

    ## scale data by train data mu/std
    test_data_mat = _utils_new._extract_adata(test_adata)
    test_adata.var['mean'] = np.mean(test_data_mat, axis=0).reshape(-1, 1)
    test_adata.var['std'] = np.std(test_data_mat, axis=0).reshape(-1, 1)
    test_data_mat = (test_data_mat - np.array(features['mean'])) / np.array(features['std'])

    y_pred = tf.nn.softmax(model.predict(test_data_mat)).numpy()
    pred_celltypes = _utils_new._prob_to_label(y_pred, encoders)
    test_adata.obs[_utils_new.PredCelltype_COLUMN] = pred_celltypes
    # change output file
    # change filepath
    test_adata.obs[['pred_celltype']].to_csv(
        "C:\\Users\\rohan\\PycharmProjects\\SciFaitEmory\\code\\predictition_files\\output.csv")

    return "done"


if __name__ == "__main__":
    print(predict())
