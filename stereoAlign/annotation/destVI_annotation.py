#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 8/28/23 4:10 PM
# @Author  : zhangchao
# @File    : destVI_annotation.py
# @Email   : zhangchao5@genomics.cn
import numpy as np
import pandas as pd
from scvi.model import CondSCVI, DestVI


def dest_vi_annotation(single_adata, spatial_adata, annotation):
    """DestVI wrapper function

    DestVI endeavors to seamlessly integrate single-cell (or single-nucleus) gene expression data with spatial gene expression data.

    :param single_adata: `AnnData`, single cell data
    :param spatial_adata: `AnnData`, spatial expression data
    :param annotation: single cell annotation cell type, which should be saved in `.obs`
    """
    assert "counts" in single_adata.layers.keys(), "can not found `counts` attribute in `.layers` keys."
    assert "counts" in spatial_adata.layers.keys(), "can not found `counts` attribute in `.layers` keys."

    # filter genes to be the same on the spatial data
    intersect = np.intersect1d(single_adata.var_names, spatial_adata.var_names)
    spatial_adata = spatial_adata[:, intersect].copy()
    single_adata = single_adata[:, intersect].copy()
    CondSCVI.setup_anndata(single_adata, layer="counts", labels_key=annotation)
    sc_model = CondSCVI(single_adata, weight_obs=False)
    sc_model.view_anndata_setup()
    sc_model.train()

    DestVI.setup_anndata(spatial_adata, layer="counts")
    st_model = DestVI.from_rna_model(spatial_adata, sc_model)
    st_model.view_anndata_setup()

    st_model.train(max_epochs=2500)

    spatial_adata.obs = pd.concat([spatial_adata.obs, st_model.get_proportions()], axis=1)

    return spatial_adata



