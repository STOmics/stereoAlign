#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 8/24/23 9:50 AM
# @Author  : zhangchao
# @File    : _desc_alignment.py
# @Email   : zhangchao5@genomics.cn
import os
import tempfile

from anndata import AnnData


def desc_alignment(adata: AnnData, batch_key, res=0.8, ncores=None, tmp_dir=None, use_gpu=False, gpu_id=None):
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
