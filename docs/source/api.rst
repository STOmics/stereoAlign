API
===


Integration
-----------

Integration method functions require the preprocessed ``anndata`` object (here ``adata``) and the name of the batch column
in ``adata.obs`` (here ``'batch'``).
The methods can be called using the following, where ``alignment`` is the name of the integration method.

For example, in order to run Harmony, on a dataset, call:

.. code-block:: python

    stereoAlign.alg.harmony_alignment(adata, batch_key="batch", n_pca=100)

Some integration methods also use cell type label as input.
For these, you need to additionally provide the corresponding label column of ``adata.obs``

.. code-block:: python

    stereoAlign.alg.scgen_alignment(adata, batch_key="batch", cell_type="celltype", epochs=200)


.. automodule:: stereoAlign.alignment

    :no-heading:

    :skip: harmony_alignment
    :skip: scanorama_alignment