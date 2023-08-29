#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 8/11/23 3:13 PM
# @Author  : zhangchao
# @File    : setup.py
# @Email   : zhangchao5@genomics.cn
import setuptools
from wheel.bdist_wheel import bdist_wheel

__version__ = "0.0.3"

with open("README.md", "r") as fh:
    long_description = fh.read()


class BDistWheel(bdist_wheel):
    def get_tag(self):
        return (self.python_tag, "none", "any")


cmdclass = {
    "bdist_wheel": BDistWheel,
}

setuptools.setup(
    name="stereoAlign",
    version=__version__,
    author="zhangchao",
    author_email="zhangchao5@genomics.cn",
    description="A toolkit package of data integration",
    long_description=long_description,
    long_description_content_type="text/markdown",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    cmdclass=cmdclass,
    python_requires='>=3.8',
    include_package_data=True,
    install_requires=[
        'matplotlib>=2.2',
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
