#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 8/24/23 9:50 AM
# @Author  : zhangchao
# @File    : _combat_alignment.py
# @Email   : zhangchao5@genomics.cn
from anndata import AnnData


def combat_alignment(adata: AnnData, batch_key):
    """ComBat wrapper function (``scanpy`` implementation)

    Using scanpy implementation of `Combat <https://scanpy.readthedocs.io/en/stable/generated/scanpy.pp.combat.html>`_

    :param adata: preprocessed ``anndata`` object
    :param batch_key: batch key in ``adata.obs``
    :return: ``anndata`` object containing the corrected feature matrix
    """
    try:
        import scanpy as sc
    except ImportError:
        raise ImportError("\nplease install scanpy:\n\n\tpip install scanpy")

    adata_int = adata.copy()
    sc.pp.combat(adata_int, key=batch_key)
    return adata_int
