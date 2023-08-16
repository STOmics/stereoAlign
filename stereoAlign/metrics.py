#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 8/16/23 2:22 PM
# @Author  : zhangchao
# @File    : metrics_utils.py
# @Email   : zhangchao5@genomics.cn
from anndata import AnnData

from stereoAlign.metrics_utils import ari, graph_connectivity, get_kbet, get_lisi, silhouette, silhouette_batch


def cal_ari(adata, cluster_key, label_key, implementation=None):
    """Adjusted Rand Index

    The adjusted rand index is a chance-adjusted rand index, which evaluates the pair-wise accuracy of clustering vs.
    ground truth label assignments.
    The score ranges between 0 and 1 with larger values indicating better conservation of the data-driven cell identity
    discovery after integration compared to annotated labels.

    :param adata: anndata object with cluster assignments in ``adata.obs[cluster_key]``
    :param cluster_key: string of column in adata.obs containing cluster assignments
    :param label_key: string of column in adata.obs containing labels
    :param implementation: if set to 'sklearn', uses sklearn's implementation,
        otherwise native implementation is taken

    This function can be applied to all integration output types.
    The ``adata`` must contain cluster assignments that are based off the knn graph given or derived from the integration
    method output.

    """
    ari(adata, cluster_key, label_key, implementation)


def cal_graph_connectivity(adata, label_key):
    """Graph Connectivity

    Quantify the connectivity of the subgraph per cell type label.
    The final score is the average for all cell type labels :math:`C`, according to the equation:

    .. math::

        GC = \\frac {1} {|C|} \\sum_{c \\in C} \\frac {|{LCC(subgraph_c)}|} {|c|}

    where :math:`|LCC(subgraph_c)|` stands for all cells in the largest connected component and :math:`|c|` stands for all cells of
    cell type :math:`c`.

    :param adata: integrated adata with computed neighborhood graph
    :param label_key: name in adata.obs containing the cell identity labels

    This function can be applied to all integration output types.
    The integrated object (``adata``) needs to have a kNN graph based on the integration output.

    """
    graph_connectivity(adata, label_key)


def cal_kbet(data, key="batch", use_rep="X_umap", alpha=0.05, n_neighbors=30):
    """K-nearest neighbour Batch Effect Test

    Calculate the K-nearest neighbors Batch Effects Test (K-BET) metric of the data regarding a specific sample attribute and embedding.
    The K-BET metric measures if cells from different samples mix well in their local neighborhood.

    :param data: Data matrix with rows for cells and columns for genes.
    :param key: The sample attribute to be consider. Must exist in ``data.obs``.
    :param use_rep: The embedding representation to be used. The key must be exist in ``data.obsm``. By default, use UMAP coordinates.
    :param alpha: Number of nearest neighbors.
    :param n_neighbors: Threshold. A cell is accepted is its K-BET p-value is greater than or equal to ``alpha``
    :return: ``stat_mean``: Mean K-BET chi-square statistic over all cells. ``pvalue_mean`` Mean K-BET  p-value over all cells. ``accept_rate`` K-BET acceptance rate of the sample.


    """

    get_kbet(data, key, use_rep, alpha, n_neighbors)


def cal_lisi(data, key="batch", use_rep="X_umap", n_neighbors=30):
    """Local inverse Simpson's Index (LISI)

    Calculate the Local inverse Simpson's Index (LISI) metric of the data regarding a specific sample attribute and embedding.
    The LISI metric measures if cells from different samples mix well in their local neighborhood.

    :param data: Data matrix with rows for cells and columns for genes.
    :param key: The sample attribute to be consider. Must exist in ``data.obs``.
    :param use_rep: The embedding representation to be used. The key must be exist in ``data.obsm``. By default, use UMAP coordinates.
    :param n_neighbors: Number of nearest neighbors.
    :return: ``lisi_mean`` Mean of calculated score. ``lower`` Lower bound of 95% confidence interval. ``upper`` Upper bound of 95% confidence interval.


    """

    get_lisi(data, key, use_rep, n_neighbors)


def cal_silhouette(adata, label_key, embed, metric="euclidean", scale=True):
    """Average silhouette width (ASW)

    Wrapper for sklearn silhouette function values range from [-1, 1] with

        * 1 indicates distinct, compact clusters
        * 0 indicates overlapping clusters
        * -1 indicates core-periphery (non-cluster) structure

    By default, the score is scaled between 0 and 1 (``scale=True``).

    :param label_key: key in adata.obs of cell labels
    :param embed: embedding key in adata.obsm, default: 'X_pca'
    :param metric: type of distance metric to use for the silhouette scores
    :param scale: default True, scale between 0 (worst) and 1 (best)

    The function requires an embedding to be stored in ``adata.obsm`` and can only be applied to feature and embedding
    integration outputs.
    Please note, that the metric cannot be used to evaluate kNN graph outputs.


    """
    silhouette(adata, label_key, embed, metric, scale)


def cal_silhouette_batch(
        adata,
        batch_key,
        label_key,
        embed,
        metric="euclidean",
        return_all=False,
        scale=True,
        verbose=True,
):
    r"""Batch ASW

    Modified average silhouette width (ASW) of batch

    This metric measures the silhouette of a given batch.
    It assumes that a silhouette width close to 0 represents perfect overlap of the batches, thus the absolute value of
    the silhouette width is used to measure how well batches are mixed.
    For all cells :math:`i` of a cell type :math:`C_j`, the batch ASW of that cell type is:

    .. math::

        batch \\, ASW_j = \\frac{1}{|C_j|} \\sum_{i \\in C_j} |silhouette(i)|

    The final score is the average of the absolute silhouette widths computed per cell type :math:`M`.

    .. math::

        batch \\, ASW = \\frac{1}{|M|} \\sum_{i \\in M} batch \\, ASW_j

    For a scaled metric (which is the default), the absolute ASW per group is subtracted from 1 before averaging, so that
    0 indicates suboptimal label representation and 1 indicates optimal label representation.

    .. math::

        batch \\, ASW_j = \\frac{1}{|C_j|} \\sum_{i \\in C_j} 1 - |silhouette(i)|

    :param batch_key: batch labels to be compared against
    :param label_key: group labels to be subset by e.g. cell type
    :param embed: name of column in adata.obsm
    :param metric: see sklearn silhouette score
    :param scale: if True, scale between 0 and 1
    :param return_all: if True, return all silhouette scores and label means
        default False: return average width silhouette (ASW)
    :param verbose: print silhouette score per group
    :return:
        Batch ASW  (always)
        Mean silhouette per group in pd.DataFrame (additionally, if return_all=True)
        Absolute silhouette scores per group label (additionally, if return_all=True)

    The function requires an embedding to be stored in ``adata.obsm`` and can only be applied to feature and embedding
    integration outputs.
    Please note, that the metric cannot be used to evaluate kNN graph outputs.


    """
    silhouette_batch(adata, batch_key, label_key, embed, metric, return_all, scale, verbose)
