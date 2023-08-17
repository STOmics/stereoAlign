#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 8/17/23 9:42 AM
# @Author  : zhangchao
# @File    : preprocessing.py
# @Email   : zhangchao5@genomics.cn
import numpy as np
import scanpy as sc
from scipy import sparse

from stereoAlign.utils import check_adata, check_batch, split_batches, merge_adata, pca_lowrank


def summarize_counts(adata, count_matrix=None, min_genes=20, min_cells=20):
    """Summarise counts of the given count matrix

    This function is useful for quality control.
    Aggregates counts per cell and per gene as well as mitochondrial fraction.
    
    Parameters
    ----------
    count_matrix: 
        count matrix, by default uses ``adata.X``
    min_cells:
        ``scanpy.pp.filter_cells`` parameter
    min_genes:
        ``scanpy.pp.filter_genes`` parameter
    
    Returns
    -------
    Include the following keys in ``adata.obs``
        'n_counts': number of counts per cell (count depth)
        'log_counts': ``np.log`` of counts per cell
        'n_genes': number of counts per gene
    """
    check_adata(adata)

    if count_matrix is None:
        count_matrix = adata.X if not sparse.issparse(adata.X) else adata.X.toarray()
    adata.obs["n_counts"] = count_matrix.sum(1)
    adata.obs["log_counts"] = np.log(adata.obs["n_counts"])
    adata.obs["n_genes"] = (count_matrix > 0).sum(1)

    sc.pp.filter_cells(adata, min_genes=min_genes)
    sc.pp.filter_genes(adata, min_cells=min_cells)


def norma_log(adata):
    """Normalization and Log transform

    :param adata:
    :return:
    """
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)


def scale_batch(adata, batch):
    """Batch-aware scaling of count matrix

    Scaling counts to a mean of 0 and standard deviation of 1 using ``scanpy.pp.scale`` for each batch separately.

    Parameters
    ----------
    adata: 
        ``anndata`` object with normalised and log-transformed counts
    batch: 
        ``adata.obs`` column

    Returns
    -------
    scaled adata
    """

    check_adata(adata)
    check_batch(batch, adata.obs)

    # Store layers for after merge (avoids vstack error in merge)
    adata_copy = adata.copy()
    tmp = dict()
    for lay in list(adata_copy.layers):
        tmp[lay] = adata_copy.layers[lay]
        del adata_copy.layers[lay]

    split = split_batches(adata_copy, batch)

    for i in split:
        sc.pp.scale(i)

    adata_scaled = merge_adata(*split, batch_key=batch, index_unique=None)
    # Reorder to original obs_name ordering
    adata_scaled = adata_scaled[adata.obs_names]

    # Add layers again
    for key in tmp:
        adata_scaled.layers[key] = tmp[key]

    del tmp
    del adata_copy

    return adata_scaled


def hvg_intersect(
        adata,
        batch,
        target_genes=2000,
        flavor="cell_ranger",
        n_bins=20,
        adataOut=False,
        n_stop=8000,
        min_genes=500,
        step_size=1000,
):
    """Highly variable gene selection

    Legacy approach to HVG selection only using HVG intersections between all batches

    Parameters
    ----------
    adata: 
        ``anndata`` object with preprocessed counts
    batch: 
        ``adata.obs`` column
    target_genes: 
        maximum number of genes (intersection reduces the number of genes)
    min_genes: 
        minimum number of intersection HVGs targeted
    step_size: 
        step size to increase HVG selection per dataset
    
    Returns 
    -------
    list of maximal ``target_genes`` number of highly variable genes
    """

    check_adata(adata)
    check_batch(batch, adata.obs)

    intersect = None
    enough = False
    n_hvg = target_genes

    split = split_batches(adata, batch)
    hvg_res = []

    for i in split:
        sc.pp.filter_genes(
            i, min_cells=1
        )  # remove genes unexpressed (otherwise hvg might break)
        hvg_res.append(
            sc.pp.highly_variable_genes(
                i, flavor=flavor, n_top_genes=n_hvg, n_bins=n_bins, inplace=False
            )
        )

    while not enough:
        genes = []

        for i in range(len(split)):
            dispersion_norm = hvg_res[i]["dispersions_norm"]
            dispersion_norm = dispersion_norm[~np.isnan(dispersion_norm)]
            dispersion_norm[::-1].sort()
            disp_cut_off = dispersion_norm[n_hvg - 1]
            gene_subset = np.nan_to_num(hvg_res[i]["dispersions_norm"]) >= disp_cut_off

            genes.append(set(split[i].var[gene_subset].index))

        intersect = genes[0].intersection(*genes[1:])
        if len(intersect) >= target_genes:
            enough = True
        else:
            if n_hvg > n_stop:
                if len(intersect) < min_genes:
                    raise Exception(
                        f"Only {len(intersect)} HVGs were found in the intersection.\n"
                        f"This is fewer than {min_genes} HVGs set as the minimum.\n"
                        "Consider raising `n_stop` or reducing `min_genes`."
                    )
                break
            n_hvg = int(n_hvg + step_size)

    if adataOut:
        return adata[:, list(intersect)].copy()

    return list(intersect)


def hvg_batch(
        adata,
        batch_key=None,
        target_genes=2000,
        flavor="cell_ranger",
        n_bins=20,
        adataOut=False,
):
    """Batch-aware highly variable gene selection

    Method to select HVGs based on mean dispersions of genes that are highly
    variable genes in all batches. Using a the top target_genes per batch by
    average normalize dispersion. If target genes still hasn't been reached,
    then HVGs in all but one batches are used to fill up. This is continued
    until HVGs in a single batch are considered.

    Parameters
    ----------
    adata:
            ``anndata`` object
    batch_key:
        ``adata.obs`` column
    target_genes:
        maximum number of genes (intersection reduces the number of genes)
    flavor:
        parameter for ``scanpy.pp.highly_variable_genes``
    n_bins:
        parameter for ``scanpy.pp.highly_variable_genes``
    adataOut:
        whether to return an ``anndata`` object or a list of highly variable genes
    """

    check_adata(adata)
    if batch_key is not None:
        check_batch(batch_key, adata.obs)

    adata_hvg = adata if adataOut else adata.copy()

    n_batches = len(adata_hvg.obs[batch_key].cat.categories)

    # Calculate double target genes per dataset
    sc.pp.highly_variable_genes(
        adata_hvg,
        flavor=flavor,
        n_top_genes=target_genes,
        n_bins=n_bins,
        batch_key=batch_key,
    )

    nbatch1_dispersions = adata_hvg.var["dispersions_norm"][
        adata_hvg.var.highly_variable_nbatches
        > len(adata_hvg.obs[batch_key].cat.categories) - 1
        ]

    nbatch1_dispersions.sort_values(ascending=False, inplace=True)

    if len(nbatch1_dispersions) > target_genes:
        hvg = nbatch1_dispersions.index[:target_genes]

    else:
        enough = False
        print(f"Using {len(nbatch1_dispersions)} HVGs from full intersect set")
        hvg = nbatch1_dispersions.index[:]
        not_n_batches = 1

        while not enough:
            target_genes_diff = target_genes - len(hvg)

            tmp_dispersions = adata_hvg.var["dispersions_norm"][
                adata_hvg.var.highly_variable_nbatches == (n_batches - not_n_batches)
                ]

            if len(tmp_dispersions) < target_genes_diff:
                print(
                    f"Using {len(tmp_dispersions)} HVGs from n_batch-{not_n_batches} set"
                )
                hvg = hvg.append(tmp_dispersions.index)
                not_n_batches += 1

            else:
                print(
                    f"Using {target_genes_diff} HVGs from n_batch-{not_n_batches} set"
                )
                tmp_dispersions.sort_values(ascending=False, inplace=True)
                hvg = hvg.append(tmp_dispersions.index[:target_genes_diff])
                enough = True

    print(f"Using {len(hvg)} HVGs")

    if not adataOut:
        del adata_hvg
        return hvg.tolist()
    else:
        return adata_hvg[:, hvg].copy()


def reduce_data(
        adata,
        pca=True,
        pca_comps=50,
        neighbors=True,
        use_rep="X_pca",
        umap=False,
):
    """Apply feature selection and dimensionality reduction steps.

    Wrapper function of PCA, neighbours computation and dimensionality
    reduction.

    Parameters
    ----------
    adata: 
        ``anndata`` object with normalised and log-transformed data in ``adata.X``
    pca: 
        whether to compute PCA
    pca_comps: 
        number of principal components
    neighbors: 
        whether to compute neighbours graph
    use_rep: 
        embedding to use for neighbourhood graph
    umap: 
        whether to compute UMAP representation
    """

    check_adata(adata)

    if pca:
        print("PCA")
        pca_lowrank(adata, n_component=pca_comps)

    if neighbors:
        print("Nearest Neigbours")
        sc.pp.neighbors(adata, use_rep=use_rep)

    if umap:
        print("UMAP")
        sc.tl.umap(adata)
