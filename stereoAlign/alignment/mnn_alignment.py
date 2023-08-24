#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 8/24/23 9:48 AM
# @Author  : zhangchao
# @File    : _mnn_alignment.py
# @Email   : zhangchao5@genomics.cn
from anndata import AnnData

from stereoAlign.utils import check_sanity, split_batches


def mnn_alignment(adata: AnnData, batch_key, hvg=None, **kwargs):
    """MNN wrapper function (``mnnpy`` implementation)

    Based on `mnnpy package <https://github.com/chriscainx/mnnpy>`_

    .. note:

        ``mnnpy`` might break with newer versions of ``numpy`` and ``pandas``

    :param adata: preprocessed ``anndata`` object
    :param batch_key: batch key in ``adata.obs``
    :param hvg: list of highly variables to subset to. If ``None``, the full dataset will be used
    :return: ``anndata`` object containing the corrected feature matrix
    """
    try:
        import mnnpy
    except ImportError:
        raise ImportError("\nplease install mnnpy:\n\n\tpip install mnnpy")

    check_sanity(adata, batch_key, hvg)
    split, categories = split_batches(adata, batch_key, return_categories=True)

    corrected, _, _ = mnnpy.mnn_correct(
        *split,
        var_subset=hvg,
        batch_key=batch_key,
        batch_categories=categories,
        index_unique=None,
        **kwargs,
    )

    return corrected
