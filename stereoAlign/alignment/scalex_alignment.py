#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 8/24/23 9:52 AM
# @Author  : zhangchao
# @File    : _scalex_alignment.py
# @Email   : zhangchao5@genomics.cn
from stereoAlign.utils import split_batches


def scalex_alignment(adata, batch_key="batch"):
    """scalex wrapper function

    Based on `scalex package <https://github.com/jsxlei/SCALEX.git>`_

    :param adata: preprocessed ``anndata`` object
    :param batch_key: batch key in ``adata.obs``
    :return: ``anndata`` object containing the corrected feature matrix
    """
    try:
        from scalex import SCALEX
    except ImportError:
        raise ImportError("\nplease install desc:\n\n\t$ pip install scalex")

    split, categories = split_batches(adata, batch_key, return_categories=True)
    corrected = SCALEX(split, processed=True, batch_name=batch_key, show=False, ignore_umap=True)

    return corrected
