#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 8/25/23 3:20 PM
# @Author  : zhangchao
# @File    : tangram_annotation.py
# @Email   : zhangchao5@genomics.cn
import scanpy as sc
import pandas as pd
import numpy as np
from typing import Union
from anndata import AnnData

from stereoAlign.utils import check_data_type


def tangram_annotation(single_adata: AnnData,
                       spatial_adata: AnnData,
                       markers: Union[None, list] = None,
                       annotation: str = "celltype",
                       mode: str = "cells",
                       device: str = "cpu",
                       perc: float = 0,
                       verbose: bool = True):
    """TANGRAM wrapper function

    Tangram endeavors to seamlessly integrate single-cell (or single-nucleus) gene expression data with spatial gene expression data.

    Based on `tangram package <https://github.com/broadinstitute/Tangram.git>`_

    The annotation results will save in `.obs`

    Parameters
    ----------
    single_adata: `AnnData`
        single cell data
    spatial_adata: `AnnData`
        spatial expression data
    markers:
        Optional. List of genes to use. If `None`, HVG genes are used.
    annotation:
        single cell annotation cell type, which should be saved in `.obs`
    mode:
        Optional. Tangram mapping mode. Currently supported: 'cell', 'clusters', 'constrained'. Default is 'cell'.
    device:
        Optional. Default is 'cpu'.
    perc:
    verbose:
    """
    try:
        import tangram as tg
    except ImportError:
        raise ImportError("\nplease install tangram:\n\n\tpip install tangram-sc")

    if markers is None:
        sc.tl.rank_genes_groups(single_adata, groupby=annotation, use_raw=False)
        marker_df = pd.DataFrame(single_adata.uns["rank_genes_groups"]["names"].iloc[0:100, :])
        gene_sc = np.unique(marker_df.melt().value.values)
        gene_st = spatial_adata.var_names.values
        markers = list(set(gene_sc).intersection(set(gene_st)))

    check_data_type(spatial_adata)

    tg.pp_adatas(single_adata, spatial_adata, genes=markers)

    assert single_adata.uns['training_genes'] == spatial_adata.uns['training_genes']

    adata_map = tg.map_cells_to_space(
        adata_sc=single_adata,
        adata_sp=spatial_adata,
        device=device,
        mode=mode,
        random_state=42,
        verbose=verbose
    )

    tg.project_cell_annotations(adata_map, spatial_adata, annotation=annotation)

    annotation_list = list(pd.unique(single_adata.obs[annotation]))
    spatial_adata.obs.drop(annotation_list, inplace=True, errors="ignore", axis=1)
    df = spatial_adata.obsm["tangram_ct_pred"][annotation_list]
    df = df.clip(df.quantile(perc), df.quantile(1 - perc), axis=1)
    df = (df - df.min()) / (df.max() - df.min())
    spatial_adata.obs = pd.concat([spatial_adata.obs, df], axis=1)
    del spatial_adata.obsm["tangram_ct_pred"]





