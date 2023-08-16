#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 8/11/23 2:53 PM
# @Author  : zhangchao
# @File    : ari.py
# @Email   : zhangchao5@genomics.cn
import numpy as np
import pandas as pd
import scipy.special
from scanpy._utils import deprecated_arg_names
from sklearn.metrics.cluster import adjusted_rand_score

from stereoAlign.utils import check_adata, check_batch


@deprecated_arg_names({"group1": "cluster_key", "group2": "label_key"})
def ari(adata, cluster_key, label_key, implementation=None):

    check_adata(adata)
    check_batch(cluster_key, adata.obs)
    check_batch(label_key, adata.obs)

    cluster_key = adata.obs[cluster_key].to_numpy()
    label_key = adata.obs[label_key].to_numpy()

    if len(cluster_key) != len(label_key):
        raise ValueError(
            f"different lengths in cluster_key ({len(cluster_key)}) and label_key ({len(label_key)})"
        )

    if implementation == "sklearn":
        return adjusted_rand_score(cluster_key, label_key)

    def binom_sum(x, k=2):
        return scipy.special.binom(x, k).sum()

    n = len(cluster_key)
    contingency = pd.crosstab(cluster_key, label_key)

    ai_sum = binom_sum(contingency.sum(axis=0))
    bi_sum = binom_sum(contingency.sum(axis=1))

    index = binom_sum(np.ravel(contingency))
    expected_index = ai_sum * bi_sum / binom_sum(n, 2)
    max_index = 0.5 * (ai_sum + bi_sum)

    return (index - expected_index) / (max_index - expected_index)
