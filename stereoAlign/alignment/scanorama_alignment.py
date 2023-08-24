#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 8/24/23 9:45 AM
# @Author  : zhangchao
# @File    : _scanorama_alignment.py
# @Email   : zhangchao5@genomics.cn
from anndata import AnnData

from stereoAlign.utils import check_sanity, split_batches, merge_adata


def scanorama_alignment(adata: AnnData, batch_key, hvg=None, **kwargs):
    """Scanorama wrapper function

    Based on `scanorama <https://github.com/brianhie/scanorama>`_

    :param adata: preprocessed ``anndata`` object
    :param batch_key: batch key in ``adata.obs``
    :param hvg: list of highly variables to subset to. If ``None``, the full dataset will be used
    :return: ``anndata`` object containing the corrected feature matrix as well as an embedding representation of the
        corrected data
    """
    try:
        import scanorama
    except ImportError:
        raise ImportError("\nplease install scanorama:\n\n\tpip install scanorama")

    check_sanity(adata, batch_key, hvg)
    split, categories = split_batches(adata.copy(), batch_key, return_categories=True)
    corrected = scanorama.correct_scanpy(split, return_dimred=True, **kwargs)
    corrected = merge_adata(
        *corrected, batch_key=batch_key, batch_categories=categories, index_unique=None
    )
    corrected.obsm["aligned_scanorama"] = corrected.obsm["X_scanorama"]

    return corrected
