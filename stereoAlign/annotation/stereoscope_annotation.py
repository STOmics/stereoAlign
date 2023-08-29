#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 8/25/23 6:07 PM
# @Author  : zhangchao
# @File    : stereoscope_annotation.py
# @Email   : zhangchao5@genomics.cn
import os
import scvi
import numpy as np
import pandas as pd
from scvi.external import RNAStereoscope, SpatialStereoscope

from stereoAlign.utils import check_data_type


def stereoscope_annotation(single_adata,
                           spatial_adata,
                           annotation,
                           train_single_model=True,
                           train_spatial_model=True,
                           model_path="./stereoscope_output"):
    """STEREOSCOPE wrapper function

    Stereoscope endeavors to seamlessly integrate single-cell (or single-nucleus) gene expression data with spatial gene expression data.

    Parameters
    ----------
    single_adata: `AnnData`
        single cell data
    spatial_adata: `AnnData`
        spatial expression data
    annotation:
        single cell annotation cell type, which should be saved in `.obs`
    train_single_model:
    train_spatial_model
    model_path:
    """
    scvi.settings.seed = 0
    check_data_type(spatial_adata)

    intersect = np.intersect1d(single_adata.var_names, spatial_adata.var_names)
    spatial_adata = spatial_adata[:, intersect].copy()
    single_adata = single_adata[:, intersect].copy()

    assert "counts" in single_adata.layers.keys(), "can not found `counts` attribute in `.layers` keys."
    RNAStereoscope.setup_anndata(single_adata, layer="counts", labels_key=annotation)

    if train_single_model:
        sc_model = RNAStereoscope(single_adata)
        sc_model.train(max_epochs=100)
        sc_model.save(os.path.join(model_path, "scmodel"), overwrite=True)
    else:
        sc_model = RNAStereoscope.load("scmodel", adata=single_adata)
        print("Loaded RNA model from file!")

    assert "counts" in spatial_adata.layers.keys(), "can not found `counts` attribute in `.layers` keys."
    SpatialStereoscope.setup_anndata(spatial_adata, layer="counts")

    if train_spatial_model:
        spatial_model = SpatialStereoscope.from_rna_model(spatial_adata, sc_model)
        spatial_model.train(max_epochs=2000)
        spatial_model.save(os.path.join(model_path, "stmodel"), overwrite=True)
    else:
        spatial_model = SpatialStereoscope.load(os.path.join(model_path, "stmodel"), adata=spatial_adata)
        print("Loaded Spatial model from file!")

    spatial_adata.obs = pd.concat([spatial_adata.obs, spatial_model.get_proportions()], axis=1)

    return spatial_adata

