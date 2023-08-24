#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 8/24/23 9:43 AM
# @Author  : zhangchao
# @File    : _harmony_alignment.py
# @Email   : zhangchao5@genomics.cn
from anndata import AnnData

from stereoAlign.utils import check_sanity, pca_lowrank


def harmony_alignment(adata: AnnData, batch_key, hvg=None, n_pca=100, **kwargs):
    """Harmony wrapper function

    Based on `harmony-pytorch <https://github.com/lilab-bcb/harmony-pytorch>`_

    :param adata: preprocessed ``anndata`` object
    :param batch_key: batch key in ``adata.obs``
    :param hvg: list of highly variables to subset to. If ``None``, the full dataset will be used
    :param n_pca: PCA component
    :return: ``anndata`` object containing the corrected feature matrix as well as an embedding representation of the corrected data
    """
    try:
        from harmony import harmonize
    except ImportError:
        raise ImportError("\nplease install harmony:\n\n\tpip install harmony-pytorch")

    check_sanity(adata, batch_key, hvg)
    pca_lowrank(adata, n_component=n_pca)
    adata.obsm["aligned_harmony"] = harmonize(adata.obsm["X_pca"], adata.obs, batch_key=batch_key)

    return adata
