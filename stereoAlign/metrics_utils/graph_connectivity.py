#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 8/11/23 2:55 PM
# @Author  : zhangchao
# @File    : graph_connectivity.py
# @Email   : zhangchao5@genomics.cn
import numpy as np
import pandas as pd
from scipy.sparse.csgraph import connected_components


def graph_connectivity(adata, label_key):

    if "neighbors" not in adata.uns:
        raise KeyError(
            "Please compute the neighborhood graph before running this function!"
        )

    clust_res = []

    for label in adata.obs[label_key].cat.categories:
        adata_sub = adata[adata.obs[label_key].isin([label])]
        _, labels = connected_components(
            adata_sub.obsp["connectivities"], connection="strong"
        )
        tab = pd.value_counts(labels)
        clust_res.append(tab.max() / sum(tab))

    return np.mean(clust_res)
