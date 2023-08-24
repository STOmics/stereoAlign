#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 8/24/23 9:49 AM
# @Author  : zhangchao
# @File    : _bbknn_alignment.py
# @Email   : zhangchao5@genomics.cn
from anndata import AnnData

from stereoAlign.utils import check_sanity, pca_lowrank


def bbknn_alignment(adata: AnnData, batch_key, hvg=None, **kwargs):
    """BBKNN wrapper function

    Based on `bbknn package <https://github.com/Teichlab/bbknn>`_

    :param adata: preprocessed ``anndata`` object
    :param batch_key: batch key in ``adata.obs``
    :param hvg: list of highly variables to subset to. If ``None``, the full dataset will be used
    :params \\**kwargs: additional parameters for BBKNN
    :return: ``anndata`` object containing the corrected graph
    """
    try:
        import bbknn
    except ImportError:
        raise ImportError("\nplease install bbknn:\n\n\tpip install bbknn")

    check_sanity(adata, batch_key, hvg)
    pca_lowrank(adata)
    if adata.n_obs < 1e5:
        return bbknn.bbknn(adata, batch_key=batch_key, copy=True, **kwargs)
    if adata.n_obs >= 1e5:
        return bbknn.bbknn(
            adata, batch_key=batch_key, neighbors_within_batch=25, copy=True, **kwargs
        )
