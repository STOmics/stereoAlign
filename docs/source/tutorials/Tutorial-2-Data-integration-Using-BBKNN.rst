Import Packages and load datasets
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: ipython3

    import os
    import sys
    
    from warnings import filterwarnings

.. code:: ipython3

    sys.path.append(os.getcwd())
    filterwarnings("ignore")

.. code:: ipython3

    import stereoAlign
    import scanpy as sc
    import numpy as np
    from anndata import AnnData

.. code:: ipython3

    dataset = AnnData.concatenate(*[sc.read_h5ad(os.path.join("demo_data", f)) for f in sorted(os.listdir("demo_data/"))])

Preprocessing datasets
~~~~~~~~~~~~~~~~~~~~~~

.. code:: ipython3

    stereoAlign.pp.summarize_counts(dataset)

.. code:: ipython3

    stereoAlign.pp.norma_log(dataset)

Optional:
~~~~~~~~~

Scaling counts to a mean of 0 and standard deviation of 1 using
``scanpy.pp.scale`` for each batch separately. It is possible to
effectively alleviate the impact of minor batch effects.

.. code:: ipython3

    dataset = stereoAlign.pp.scale_batch(dataset, batch="batch")

Visualization original dataset
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Before proceeding with the integration process, it is essential to
visualize the original dataset in its unadulterated form. This
visualization serves as a valuable reference point for understanding the
inherent characteristics and structure of the data.

.. code:: ipython3

    stereoAlign.pp.reduce_data(dataset, pca=True, pca_comps=100, neighbors=True, use_rep="X_pca", umap=True)


.. parsed-literal::

    PCA
    Nearest Neigbours
    UMAP


.. code:: ipython3

    sc.tl.leiden(dataset)

.. code:: ipython3

    sc.pl.umap(dataset, color=["batch", "leiden", "celltype"])



.. image:: Tutorial-2-Data-integration-Using-BBKNN_files/Tutorial-2-Data-integration-Using-BBKNN_13_0.png


Calculation metric score using original datasets
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: ipython3

    stat_mean, pvalue_mean, accept_rate = stereoAlign.metrics.kbet(
        dataset, 
        key="batch", 
        use_rep="X_umap", 
        alpha=0.1,
        n_neighbors=15)

.. code:: ipython3

    stat_mean, pvalue_mean, accept_rate




.. parsed-literal::

    (17.377788396500012, 0.08517801527945788, 0.19493427380570696)



.. code:: ipython3

    stereoAlign.metrics.graph_connectivity(dataset, label_key="celltype")




.. parsed-literal::

    0.9917491749174918



.. code:: ipython3

    stereoAlign.metrics.silhouette(dataset, label_key="celltype", embed="X_umap")




.. parsed-literal::

    0.5771105885505676



Using ``BBKNN`` to integration datasets
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: ipython3

    bbknn_dataset = dataset.copy()
    bbknn_corrected = stereoAlign.alg.bbknn_alignment(bbknn_dataset, batch_key="batch")


.. parsed-literal::

    /home/zhangchao/anaconda3/envs/py38/lib/python3.8/site-packages/faiss/loader.py:28: DeprecationWarning: distutils Version classes are deprecated. Use packaging.version instead.
      if LooseVersion(numpy.__version__) >= "1.19":
    /home/zhangchao/anaconda3/envs/py38/lib/python3.8/site-packages/setuptools/_distutils/version.py:345: DeprecationWarning: distutils Version classes are deprecated. Use packaging.version instead.
      other = LooseVersion(other)


.. code:: ipython3

    bbknn_corrected




.. parsed-literal::

    AnnData object with n_obs × n_vars = 3119 × 20341
        obs: 'celltype', 'batch', 'n_counts', 'log_counts', 'n_genes', 'leiden'
        var: 'n_cells', 'mean-0', 'std-0', 'mean-1', 'std-1', 'mean-2', 'std-2'
        uns: 'pca', 'neighbors', 'umap', 'leiden', 'batch_colors', 'leiden_colors', 'celltype_colors'
        obsm: 'spatial', 'X_pca', 'X_umap', 'X_umap_knn_connectivity', 'X_umap_knn_distances'
        obsp: 'distances', 'connectivities'



Visualization of integrated datasets
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: ipython3

    stereoAlign.pp.reduce_data(
        bbknn_corrected, pca=True, pca_comps=100, neighbors=False, use_rep="X_pca", umap=True)


.. parsed-literal::

    PCA
    UMAP


.. code:: ipython3

    sc.tl.leiden(bbknn_corrected)

.. code:: ipython3

    sc.pl.umap(bbknn_corrected, color=["batch", "leiden", "celltype"])



.. image:: Tutorial-2-Data-integration-Using-BBKNN_files/Tutorial-2-Data-integration-Using-BBKNN_25_0.png


Calculation metirc score using integrated datasets by ``BBKNN``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

| In order to quantify the performance and efficacy of the integrated
  datasets using the powerful ``BBKNN``, we can calculate the metric
  score.
| This score serves as a reliable measure of the harmonious blending and
  seamless integration achieved through the amalgamation process.

-  Calculate the K-nearest neighbors Batch Effects Test (K-BET) metric
   of the data regarding a specific sample attribute and embedding. The
   K-BET metric measures if cells from different samples mix well in
   their local neighborhood.

-  Should the p-value surpass the predetermined alpha threshold, it is
   indicative that the data batch effect has indeed been successfully
   eradicated. This serves as a testament to the efficacy of the applied
   methodologies and techniques employed in the removal process.

.. code:: ipython3

    stat_mean, pvalue_mean, accept_rate = stereoAlign.metrics.kbet(
        bbknn_corrected, 
        key="batch", 
        use_rep="X_umap", 
        alpha=0.05,
        n_neighbors=15)

.. code:: ipython3

    stat_mean, pvalue_mean, accept_rate




.. parsed-literal::

    (7.685325943884262, 0.18556349233809058, 0.5338249438922732)



-  Calculate the Local inverse Simpson's Index (LISI) metric of the data
   regarding a specific sample attribute and embedding. The LISI metric
   measures if cells from different samples mix well in their local
   neighborhood.

-  The larger the ``ilisi_mean`` value, the better.

.. code:: ipython3

    ilisi_mean, lower, upper = stereoAlign.metrics.lisi(
        bbknn_corrected, 
        key="batch", 
        use_rep="X_umap", 
        n_neighbors=15)

.. code:: ipython3

    ilisi_mean, lower, upper




.. parsed-literal::

    (2.099519537611782, 2.0818878012324475, 2.1171512739911167)



-  The smaller the ``clisi_mean`` value, the better.

.. code:: ipython3

    clisi_mean, lower, upper = stereoAlign.metrics.lisi(
        bbknn_corrected, 
        key="celltype", 
        use_rep="X_umap", 
        n_neighbors=15)

.. code:: ipython3

    clisi_mean, lower, upper




.. parsed-literal::

    (1.2433578961733485, 1.2310056237428135, 1.2557101686038834)



-  Quantify the connectivity of the subgraph per cell type label.

.. code:: ipython3

    stereoAlign.metrics.graph_connectivity(bbknn_corrected, label_key="celltype")




.. parsed-literal::

    0.9910141958232087



-  Average silhouette width (ASW)

   The values range from [-1, 1] with

   ::

       * 1 indicates distinct, compact clusters
       * 0 indicates overlapping clusters
       * -1 indicates core-periphery (non-cluster) structure

   By default, the score is scaled between 0 and 1 (``scale=True``).

.. code:: ipython3

    stereoAlign.metrics.silhouette(bbknn_corrected, label_key="celltype", embed="X_umap")




.. parsed-literal::

    0.6331341564655304



-  Modified average silhouette width (ASW) of batch

   This metric measures the silhouette of a given batch. It assumes that
   a silhouette width close to 0 represents perfect overlap of the
   batches, thus the absolute value of the silhouette width is used to
   measure how well batches are mixed.

.. code:: ipython3

    stereoAlign.metrics.silhouette_batch(
        bbknn_corrected, batch_key="batch", label_key="celltype", embed="X_umap")


.. parsed-literal::

    mean silhouette per group:                                 silhouette_score
    group                                           
    glomerular layer (GL)                   0.721066
    granule cell layer (GCL)                0.872616
    mitral cell layer (MCL)                 0.845666
    olfactory nerve layer (ONL)             0.872244
    rostral migratory stream (RMS)          0.698698




.. parsed-literal::

    0.8020580056902566



