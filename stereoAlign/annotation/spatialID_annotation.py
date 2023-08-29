#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 8/28/23 4:09 PM
# @Author  : zhangchao
# @File    : spatialid_annotation.py
# @Email   : zhangchao5@genomics.cn
import os.path


def spatial_id_annotation(spatial_data,
                          single_data=None,
                          annotation=None,
                          markers=None,
                          pca_dim=200,
                          k_graph=15,
                          w_cls=20,
                          w_dae=1,
                          w_gae=1,
                          gpu="0",
                          dnn_pt=None,
                          model_path="./output",):
    try:
        import spatialID
    except ImportError:
        raise ImportError("\nplease install spatialID:\n\n\tpip install spatialID")


