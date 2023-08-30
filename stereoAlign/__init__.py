#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 8/11/23 2:06 PM
# @Author  : zhangchao
# @File    : __init__.py.py
# @Email   : zhangchao5@genomics.cn
from . import alignment, metrics, preprocessing, annotation
from warnings import filterwarnings

filterwarnings("ignore")
__version__ = "0.0.2"

alg = alignment
ann = annotation
pp = preprocessing
