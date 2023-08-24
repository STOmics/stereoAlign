#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 8/24/23 9:46 AM
# @Author  : zhangchao
# @File    : _scgen_alignment.py
# @Email   : zhangchao5@genomics.cn
from anndata import AnnData

from stereoAlign.utils import check_sanity


def scgen_alignment(adata: AnnData, batch_key, cell_type, epochs=200, hvg=None, **kwargs):
    """scGen wrapper function

    Based on `scgen`_ with parametrization taken from the tutorial `notebook`_.

    .. _scgen: https://github.com/theislab/scgen
    .. _notebook: https://scgen.readthedocs.io/en/stable/tutorials/scgen_batch_removal.html

    :param adata: preprocessed ``anndata`` object
    :param batch_key: batch key in ``adata.obs``
    :param hvg: list of highly variables to subset to. If ``None``, the full dataset will be used
    :return: ``anndata`` object containing the corrected feature matrix
    """
    try:
        from scgen import SCGEN
    except ImportError:
        raise ImportError("\nplease install scgen:\n\n\tpip install scgen")

    check_sanity(adata, batch_key, hvg)

    net_adata = adata.copy()
    if hvg is not None:
        net_adata = net_adata[:, hvg].copy()

    SCGEN.setup_anndata(net_adata, batch_key=batch_key, labels_key=cell_type)
    model = SCGEN(net_adata)
    model.train(
        max_epochs=epochs,
        batch_size=32,
        early_stopping=True,
        early_stopping_patience=25,
    )
    corrected_adata = model.batch_removal(**kwargs)
    return corrected_adata
