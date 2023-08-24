#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 8/24/23 9:47 AM
# @Author  : zhangchao
# @File    : _scvi_alignment.py
# @Email   : zhangchao5@genomics.cn
from anndata import AnnData

from stereoAlign.utils import check_sanity


def scvi_alignment(adata: AnnData, batch_key, hvg=None, return_model=False, max_epochs=None):
    """scVI wrapper function

    Based on scvi-tools version >=0.16.0 (available through `conda <https://docs.scvi-tools.org/en/stable/installation.html>`_)

    .. note::
        scVI expects only non-normalized (count) data on highly variable genes!

    :param adata: preprocessed ``anndata`` object
    :param batch_key: batch key in ``adata.obs``
    :param hvg: list of highly variables to subset to. If ``None``, the full dataset will be used
    :return: ``anndata`` object containing the corrected feature matrix as well as an embedding representation of the
        corrected data
    """
    try:
        from scvi.model import SCVI
    except ImportError:
        raise ImportError("\nplease install scvi:\n\n\tpip install scvi-tools")

    check_sanity(adata, batch_key, hvg)

    # Check for counts data layer
    if "counts" not in adata.layers:
        raise TypeError(
            "Adata does not contain a `counts` layer in `adata.layers[`counts`]`"
        )

    # Defaults from SCVI github tutorials scanpy_pbmc3k and harmonization
    n_latent = 30
    n_hidden = 128
    n_layers = 2

    # copying to not return values added to adata during setup_anndata
    net_adata = adata.copy()
    if hvg is not None:
        net_adata = adata[:, hvg].copy()
    SCVI.setup_anndata(net_adata, layer="counts", batch_key=batch_key)

    vae = SCVI(
        net_adata,
        gene_likelihood="nb",
        n_layers=n_layers,
        n_latent=n_latent,
        n_hidden=n_hidden,
    )
    train_kwargs = {"train_size": 1.0}
    if max_epochs is not None:
        train_kwargs["max_epochs"] = max_epochs
    vae.train(**train_kwargs)
    adata.obsm["aligned_scvi"] = vae.get_latent_representation()

    if not return_model:
        return adata
    else:
        return vae