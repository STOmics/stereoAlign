#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 8/24/23 9:51 AM
# @Author  : zhangchao
# @File    : _spatialign_alignment.py
# @Email   : zhangchao5@genomics.cn
from anndata import AnnData


def spatialign_alignment(adata: AnnData, batch_key="batch", latent_dims=100, n_neigh=15, is_verbose=False, tau1=0.2,
                         tau2=1., tau3=0.5, save_path="./spatialign_output"):
    """spatiAlign wrapper function

     Based on `sptiAlign package <https://github.com/zhangchao162/Spatialign.git>`_

    :param adata: preprocessed ``anndata`` object
    :param batch_key: batch key in ``adata.obs``
    :param latent_dims: The number of embedding dimensions, default, 100.
    :param n_neigh: The number of neighbors selected when constructing a spatial neighbor graph. default, 15.
    :param is_verbose: Whether the detail information is print, default, True.
    :param tau1: Instance level and pseudo prototypical cluster level contrastive learning parameters, default, 0.2
    :param tau2: Pseudo prototypical cluster entropy parameter, default, 1.
    :param tau3: Cross-batch instance self-supervised learning parameter, default, 0.5
    :param save_path: The path of alignment dataset and saved spatialign.
    :return: ``anndata`` object containing the corrected feature matrix
    """

    try:
        from spatialign import Spatialign
    except ImportError:
        raise ImportError("\nplease install desc:\n\n\t$ pip install spatialign==0.0.2a0")

    model = Spatialign(
        merge_data=adata, batch_key=batch_key, is_hvg=False, is_reduce=False, n_pcs=100, n_hvg=2000, n_neigh=n_neigh,
        is_undirected=True, latent_dims=latent_dims, is_verbose=is_verbose, seed=42, gpu=None, save_path=save_path
    )
    model.train(lr=1e-3, max_epoch=500, alpha=0.5, patient=15, tau1=tau1, tau2=tau2, tau3=tau3)

    corrected_data = model.alignment()

    return corrected_data
