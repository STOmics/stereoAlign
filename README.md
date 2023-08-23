[![python >3.8.8](https://img.shields.io/badge/python-3.8.8-brightgreen)](https://www.python.org/)
[![Downloads](https://static.pepy.tech/badge/stereoalign)](https://pepy.tech/project/stereoalign)
[![Downloads](https://static.pepy.tech/badge/stereoalign/month)](https://pepy.tech/project/stereoalign)
[![Downloads](https://static.pepy.tech/badge/stereoalign/week)](https://pepy.tech/project/stereoalign)
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
The `stereoAlign` python package is available on [Pypi](https://pypi.org/project/stereoAlign/) and can be installed through      
```python
pip install stereoAlign
```            
Import `stereoAlign` in python  
```python
import stereoAlign
```
We created the python package called `stereoAlign` that uses `scanpy` to streamline the integration of single-cell datasets and evaluate the results. The package contains several modules for preprocessing an `anndata` object, running integration methods and evaluating the resulting using a number of metrics. Functions for the data integration methods are in `stereoAlign.alignment` or for short `stereoAlign.alg` and metrics are under `stereoAlign.metrics`.
       
# Tutorials        
Quick Start [https://stereoalign-tutorial.readthedocs.io/en/latest/index.html#](https://stereoalign-tutorial.readthedocs.io/en/latest/index.html#)        
        
# Integration Tools        
This repository contains the code for the `stereoAlign` package for data integration tools.           
This toolkit that is compared include:          

| **Integration Method** |               **Invoke**               |      **Output**      |                                     **Recommendation index**                                      |                                                                          **Describe**                                                                           |  Language  |
|:----------------------:|:--------------------------------------:|:--------------------:|:-------------------------------------------------------------------------------------------------:|:---------------------------------------------------------------------------------------------------------------------------------------------------------------:|:----------:|
|       `Harmony`        |  `stereoAlign.alg.harmony_alignment`   |      Embedding       |              <font color=YellowGreen size=5 face=Rockwell Extra Bold>__***__</font>               |                                   Leverage iterative clustering with maximum diversity for batch correction and integration.                                    | Python / R |
|      `Scanorama`       | `stereoAlign.alg.scanorama_alignment`  |      Embedding       |              <font color=YellowGreen size=5 face=Rockwell Extra Bold>__****__</font>              |                            Leverage computer vision algorithms for panorama stitching matches for batch correction and integration.                             |   Python   |
|        `scGEN`         |   `stereoAlign.alg.scgen_alignment`    |       Features       |              <font color=YellowGreen size=5 face=Rockwell Extra Bold>__****__</font>              |           Use variational autoencoders to reduce the dimension of gene expression matrix, and apply extra label to batch correction and integrations.           |   Python   |
|         `scvi`         |    `stereoAlign.alg.scvi_alignment`    |      Embedding       |              <font color=YellowGreen size=5 face=Rockwell Extra Bold>__****__</font>              |                                                  Use variational autoencoders to batch correct and integrate.                                                   |   Python   |
|         `MNN`          |    `stereoAlign.alg.mnn_alignment`     |       Features       |               <font color=YellowGreen size=5 face=Rockwell Extra Bold>__**__</font>               |                                             Leverage mutual nearest neighbor for batch correction and integration.                                              |   Python   |                                                                                                                                                                                                                  
|        `BBKNN`         |   `stereoAlign.alg.bbknn_alignment`    |        Graph         |              <font color=YellowGreen size=5 face=Rockwell Extra Bold>__***__</font>               |     BBKNN taking each cell and identifying a (smaller) k nearest neighbours in each batch separately, and merged into a final neighbour list for the cell.      |   Python   |                                                                                                                                                                                                                  
|        `Combat`        |   `stereoAlign.alg.combat_alignment`   |       Features       |               <font color=YellowGreen size=5 face=Rockwell Extra Bold>__*__</font>                |                                             Leverage empirical Bayesian model for batch correction and integration.                                             |   Python   |                                                                                                                                                                                                                 
|         `DESC`         |    `stereoAlign.alg.desc_alignment`    |      Embedding       |               <font color=YellowGreen size=5 face=Rockwell Extra Bold>__**__</font>               | DESC eventually eliminates batch effects by recurrent self-learning, as long as technological changes between batches are fewer than real biological variances. |   Python   |                                                                                                                                                                                                                  
|       `PRECAST`        |                                        |      Embedding       |              <font color=YellowGreen size=5 face=Rockwell Extra Bold>__****__</font>              |                            PRECAST is a probabilistic model-based approach that integrates SRT datasets from multiple tissue slides.                            |     R      |                                                                                                                                                                                                                  
|      `spatiAlign`      | `stereoAlign.alg.spatialign_alignment` | Embedding / Features |             <font color=YellowGreen size=5 face=Rockwell Extra Bold>__*****__</font>              |                       spatiAlign is an unsupervised across domain adapatation methods for SRT datasets batch correction and integration.                        |   Python   |
|        `SCALEX`        |   `stereoAlign.alg.scalex_alignment`   |      Embedding       |               <font color=YellowGreen size=5 face=Rockwell Extra Bold>__**__</font>               |                                    SCALEX leverage domain specific batch normalization for batch correction and integration.                                    |   Python   |


        
# Metric        
|  **Metric Method**  |                  **Invoke**                   |
|:-------------------:|:---------------------------------------------:|
|         ARI         |         `stereoAlign.metrics.cal_ari`         |
| Graph connectivity  | `stereoAlign.metrics.cal_graph_connectivity`  |
|        KBET         |        `stereoAlign.metrics.cal_kbet`         |
|        LISI         |        `stereoAlign.metrics.cal_lisi`         |
|     Silhouette      |     `stereoAlign.metrics.cal_silhouette`      |
        

 
# Disclaimer        
***This is not an official product.***    
            
            