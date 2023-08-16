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
    r"""Graph Connectivity

    Quantify the connectivity of the subgraph per cell type label.
    The final score is the average for all cell type labels :math:`C`, according to the equation:

    .. math::

        GC = \\frac {1} {|C|} \\sum_{c \\in C} \\frac {|{LCC(subgraph_c)}|} {|c|}

    where :math:`|LCC(subgraph_c)|` stands for all cells in the largest connected component and :math:`|c|` stands for all cells of
    cell type :math:`c`.

    Parameters
    ----------
    adata:
        integrated adata with computed neighborhood graph
    label_key:
        name in adata.obs containing the cell identity labels

    This function can be applied to all integration output types.
    The integrated object (``adata``) needs to have a kNN graph based on the integration output.

    """
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
