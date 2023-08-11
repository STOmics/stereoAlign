#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 8/11/23 3:13 PM
# @Author  : zhangchao
# @File    : setup.py
# @Email   : zhangchao5@genomics.cn
import setuptools

__version__ = "0.0.1"

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="stereoAlign",
    version=__version__,
    author="zhangchao",
    author_email="zhangchao5@genomics.cn",
    description="A toolkit package of batch effects removal",
    long_description=long_description,
    long_description_content_type="text/markdown",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.8',
    include_package_data=True,
    install_requires=[
        'matplotlib>=2.2'
        'tensorflow',
        'scanpy',
        'louvain',
        'python-igraph',
        'h5py',
        'pandas',
        'numpy',
        'pandas',
        'scanpy',
        'scikit-learn',
        'harmony-pytorch',
        'scanorama',
        'scgen',
        'scvi-tools',
        'mnnpy',
        'bbknn',
        'desc'
    ],
)