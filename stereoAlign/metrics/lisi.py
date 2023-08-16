#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 8/11/23 2:59 PM
# @Author  : zhangchao
# @File    : lisi.py
# @Email   : zhangchao5@genomics.cn
import psutil
import pandas as pd
import numpy as np
from anndata import AnnData
from functools import partial
from multiprocessing import Pool

from stereoAlign.utils import get_neighbors


def one_sample_lisi(neighbors_indices, attr_values):
    probability = pd.Series(attr_values[neighbors_indices]).value_counts(normalize=True, sort=False).values
    score = 1 / (probability ** 2).sum()
    return score


def get_lisi(
        data: AnnData,
        key: str = "batch",
        use_rep: str = "X_umap",
        n_neighbors: int = 30):
    """Local inverse Simpson's Index (LISI)

    Calculate the Local inverse Simpson's Index (LISI) metric of the data regarding a specific sample attribute and embedding.
    The LISI metric measures if cells from different samples mix well in their local neighborhood.

    Parameters
    ----------
    data:
        Data matrix with rows for cells and columns for genes.
    key:
        The sample attribute to be consider. Must exist in ``data.obs``.
    use_rep:
        The embedding representation to be used. The key must be exist in ``data.obsm``. By default, use UMAP coordinates.
    n_neighbors:
        Number of nearest neighbors.

    Returns
    -------
    lisi_mean:
        Mean of calculated score.
    lower:
        Lower bound of 95% confidence interval.
    upper:
        Upper bound of 95% confidence interval.

    """

    assert key in data.obs_keys()

    n_sample = data.shape[0]
    get_neighbors(data, n_neighbors=n_neighbors, use_rep=use_rep)

    # add itself into the knn connectivity graph
    assert f"{use_rep}_knn_connectivity" in data.obsm_keys(), f"Error, can not found '{use_rep}_knn_connectivity' in .obsm_keys(). Please calculate nearest neighbors graph first."
    indices = np.concatenate((np.arange(n_sample).reshape(-1, 1), data.obsm[f"{use_rep}_knn_connectivity"][:, :-1]),
                             axis=1)

    partial_lisi = partial(
        one_sample_lisi,
        attr_values=data.obs[key].values.copy())

    n_jobs = psutil.cpu_count(logical=False)
    if n_jobs is None:
        n_jobs = psutil.cpu_count(logical=True)

    with Pool(n_jobs) as p:
        results = p.map(partial_lisi, indices)
    results = np.array(results)
    lisi_mean = results.mean()

    std = results.std()
    lower = lisi_mean - 1.96 * std / np.sqrt(n_sample)
    upper = lisi_mean + 1.96 * std / np.sqrt(n_sample)
    return lisi_mean, lower, upper
