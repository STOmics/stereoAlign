Installation
==============

We provide two ways to install the package of stereoAlign.

Please note that the current stereoAlign version offers full support of Linux operating system. Further version for other operating systems would be released soon.


1. Python
----------------
Install through `Pypi <https://pypi.org/project/stereoAlign/>`_.

.. code-block:: python

	pip install stereoAlign

or

.. code-block:: python

    git clone https://github.com/STOmics/stereoAlign.git

    cd stereoAlign

    python setup.py install


2. Anaconda
---------------
For convenience, we suggest using a separate conda environment for running stereoAlign. Please ensure Annaconda3 is installed.

Create conda environment and install stereoAlign package.

.. code-block:: python

   # create an environment called Spatialign
   conda create -n stereoAlign python=3.8

   # activate your environment
   conda activate stereoAlign

   # install package

   pip install stereoAlign

   or

   git clone https://github.com/STOmics/stereoAlign.git

   cd stereoAlign

   python setup.py build

   python setup.py install --user

   #To use the environment in jupyter notebook, add python kernel for this environment.

   pip install ipykernel

   python -m ipykernel install --user --name=stereoAlign