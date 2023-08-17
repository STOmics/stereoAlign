#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 8/11/23 2:08 PM
# @Author  : zhangchao
# @File    : alignment.py
# @Email   : zhangchao5@genomics.cn
import os
import tempfile
from anndata import AnnData

from stereoAlign.utils import check_sanity, pca_lowrank, split_batches, merge_adata


def harmony_alignment(adata: AnnData, batch_key: str, hvg=None, n_pca: int = 100, **kwargs):
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


def scanorama_alignment(adata: AnnData, batch_key: str, hvg=None, **kwargs):
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


def scgen_alignment(adata: AnnData, batch_key: str, cell_type: str, epochs: int = 200, hvg=None, **kwargs):
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


def scvi_alignment(adata: AnnData, batch_key: str, hvg=None, return_model=False, max_epochs=None):
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


def mnn_alignment(adata: AnnData, batch_key: str, hvg=None, **kwargs):
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


def bbknn_alignment(adata: AnnData, batch_key: str, hvg=None, **kwargs):
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


def combat_alignment(adata: AnnData, batch_key: str):
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


def desc_alignment(adata: AnnData, batch_key: str, res=0.8, ncores=None, tmp_dir=None, use_gpu=False, gpu_id=None):
    """DESC wrapper function

    Based on `desc package <https://github.com/eleozzr/desc>`_
    Parametrization was taken from: https://github.com/eleozzr/desc/issues/28 as suggested by the developer (rather
    than from the tutorial notebook).

    :param adata: preprocessed ``anndata`` object
    :param batch_key: batch key in ``adata.obs``
    :return: ``anndata`` object containing the corrected embedding
    """
    try:
        import desc
    except ImportError:
        raise ImportError("\nplease install desc:\n\n\t$ pip install desc")

    if tmp_dir is None:
        temp_dir = tempfile.TemporaryDirectory()
        tmp_dir = temp_dir.name

    # Set number of CPUs to all available
    if ncores is None:
        ncores = os.cpu_count()

    adata_out = adata.copy()

    adata_out = desc.scale_bygroup(adata_out, groupby=batch_key, max_value=6)

    adata_out = desc.train(
        adata_out,
        dims=[adata.shape[1], 128, 32],
        tol=0.001,
        n_neighbors=10,
        batch_size=256,
        louvain_resolution=res,
        save_encoder_weights=False,
        save_dir=tmp_dir,
        do_tsne=False,
        use_GPU=use_gpu,
        GPU_id=gpu_id,
        num_Cores=ncores,
        use_ae_weights=False,
        do_umap=False,
    )

    adata_out.obsm["aligned_desc"] = adata_out.obsm["X_Embeded_z" + str(res)]

    return adata_out
