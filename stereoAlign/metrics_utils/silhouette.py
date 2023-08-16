#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 8/11/23 3:00 PM
# @Author  : zhangchao
# @File    : silhouette.py
# @Email   : zhangchao5@genomics.cn
import numpy as np
import pandas as pd
from scanpy._utils import deprecated_arg_names
from sklearn.metrics.cluster import silhouette_samples, silhouette_score


@deprecated_arg_names({"group_key": "label_key"})
def silhouette(adata, label_key, embed, metric="euclidean", scale=True):

    if embed not in adata.obsm.keys():
        print(adata.obsm.keys())
        raise KeyError(f"{embed} not in obsm")
    asw = silhouette_score(
        X=adata.obsm[embed], labels=adata.obs[label_key], metric=metric
    )
    if scale:
        asw = (asw + 1) / 2
    return asw


@deprecated_arg_names({"group_key": "label_key"})
def silhouette_batch(
        adata,
        batch_key,
        label_key,
        embed,
        metric="euclidean",
        return_all=False,
        scale=True,
        verbose=True,
):

    if embed not in adata.obsm.keys():
        print(adata.obsm.keys())
        raise KeyError(f"{embed} not in obsm")

    sil_per_label = []
    for group in adata.obs[label_key].unique():
        adata_group = adata[adata.obs[label_key] == group]
        n_batches = adata_group.obs[batch_key].nunique()

        if (n_batches == 1) or (n_batches == adata_group.shape[0]):
            continue

        sil = silhouette_samples(
            adata_group.obsm[embed], adata_group.obs[batch_key], metric=metric
        )

        # take only absolute value
        sil = [abs(i) for i in sil]

        if scale:
            # scale s.t. highest number is optimal
            sil = [1 - i for i in sil]

        sil_per_label.extend([(group, score) for score in sil])

    sil_df = pd.DataFrame.from_records(
        sil_per_label, columns=["group", "silhouette_score"]
    )

    if len(sil_per_label) == 0:
        sil_means = np.nan
        asw = np.nan
    else:
        sil_means = sil_df.groupby("group").mean()
        asw = sil_means["silhouette_score"].mean()

    if verbose:
        print(f"mean silhouette per group: {sil_means}")

    if return_all:
        return asw, sil_means, sil_df

    return asw
