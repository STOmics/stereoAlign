[![python >3.8.8](https://img.shields.io/badge/python-3.8.8-brightgreen)](https://www.python.org/)          
# A toolkit for data integration
***        
Behold, a magnificent toolkit emerges, designed to seamlessly integrate the vast realm of data in the realm of single-cell genomics.   Within its hallowed digital confines, this toolkit bestows upon researchers the power to harmoniously merge and analyze the intricate tapestry of cellular information. Crafted with meticulous care, this toolkit embodies the epitome of elegance and efficiency.   It serves as a conduit, enabling the synthesis of diverse datasets from single-cell genomics experiments, unlocking new realms of knowledge and understanding. The power of data integration within this toolkit extends beyond mere aggregation.   It empowers researchers to unravel the intricate web of cellular interactions, uncovering hidden patterns, identifying novel cell types, and discerning the complex dynamics that govern cellular behavior. Through its refined algorithms and advanced statistical techniques, this toolkit illuminates the path towards deeper insights and discoveries.   It refines and enhances the quality of data, mitigating confounding factors and removing noise, thus unveiling the true essence of the cellular landscape. As the sun sets on the horizon of single-cell genomics, this toolkit emerges as a guiding beacon, illuminating the path towards a more comprehensive understanding of cellular complexity.   Embrace its power and unlock the secrets hidden within the realm of single-cell genomics.          
            
# Dependence        
[![torch-1.10.0](https://img.shields.io/badge/torch-1.10.0-red)](https://pytorch.org/get-started/previous-versions/)
[![pandas-1.2.4](https://img.shields.io/badge/pandas-1.2.4-lightgrey)](https://github.com/pandas-dev/pandas)
[![scikit-learn-0.24](https://img.shields.io/badge/scikit-0.24.x-brightgreen)](https://github.com/scikit-learn/scikit-learn/tree/0.24.X)
[![scipy-0.12.x](https://img.shields.io/badge/scipy-0.12.x-yellow)](https://github.com/scipy/scipy/tree/maintenance/0.12.x)
[![scanpy-1.9.1](https://img.shields.io/badge/scanpy-1.9.1-informational)](https://pypi.org/project/scanpy/)           
                
# Install         
The `stereoAlign` python package is available on [Pypi]() and can be installed through      
```python
pip install stereoAlign
```            
Import `stereoAlign` in python  
```python
import stereoAlign
```
We created the python package called `stereoAlign` that uses `scanpy` to streamline the integration of single-cell datasets and evaluate the results. The package contains several modules for preprocessing an `anndata` object, running integration methods and evaluating the resulting using a number of metrics. Functions for the data integration methods are in `stereoAlign.alignment` or for short `stereoAlign.alg` and metrics are under `stereoAlign.metrics`.
       
        
# Integration Tools        
This repository contains the code for the `stereoAlign` package for data integration tools.           
This toolkit that is compared include:          

| **Integration Method** |                  **Invoke**                  |                                                                                                   **Describe**                                                                                                   |
|:----------------------:|:--------------------------------------------:|:----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------:|
|        Harmony         |     `stereoAlign.alg.harmony_alignment`      |                                                        Harmony wrapper function based on `harmony-pytorch` <https://github.com/lilab-bcb/harmony-pytorch>                                                        |
|       Scanorama        |    `stereoAlign.alg.scanorama_alignment`     |                                                             Scanorama wrapper function based on `scanorama` <https://github.com/brianhie/scanorama>                                                              |
|         scGEN          |      `stereoAlign.alg.scgen_alignment`       | scGen wrapper function based on `scgen` <https://github.com/theislab/scgen> with parametrization taken from the tutorial `notebook` <https://scgen.readthedocs.io/en/stable/tutorials/scgen_batch_removal.html>. |
|          scvi          |       `stereoAlign.alg.scvi_alignment`       |                                     scVI wrapper function based on scvi-tools version >=0.16.0 (available through <https://docs.scvi-tools.org/en/stable/installation.html>)                                     |
|          MNN           |       `stereoAlign.alg.mnn_alignment`        |                                                  MNN wrapper function (``mnnpy`` implementation) based on `mnnpy package` <https://github.com/chriscainx/mnnpy>                                                  |
|         BBKNN          |      `stereoAlign.alg.bbknn_alignment`       |                                                               BBKNN wrapper function based on `bbknn package` <https://github.com/Teichlab/bbknn>                                                                |
|         Combat         |      `stereoAlign.alg.combat_alignment`      |                      ComBat wrapper function (``scanpy`` implementation) using scanpy implementation of `Combat` <https://scanpy.readthedocs.io/en/stable/generated/scanpy.pp.combat.html>                       |
|          DESC          |       `stereoAlign.alg.desc_alignment`       |                                                                 DESC wrapper function based on `desc package` <https://github.com/eleozzr/desc>                                                                  |
|      More methods      | Additional methods are being incorporated... |                                                                                                                                                                                                                  |

        
# Metric        
|  **Metric Method**  | **Invoke**                               |
|:-------------------:|:-----------------------------------------:|
|         ARI         | `stereoAlign.metrics.ari`                |
| Graph connectivity  | `stereoAlign.metrics.graph_connectivity` |
|        KBET         | `stereoAlign.metrics.get_kbet`           |
|        LISI         | `stereoAlign.metrics.get_lisi`           |
|     Silhouette      | `stereoAlign.metrics.silhouette`         |
        

 
# Disclaimer        
***This is not an official product.***    
            
            
