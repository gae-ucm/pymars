
pyMARS
======

pyMARS is a python-based package under active development for IACT event analysis of the MAGIC telescopes. The purpose of this repo is to quickly and easily visualize and process the output of `CTLearn <https://github.com/ctlearn-project/ctlearn>`_ using data from the MAGIC telescopes. pyMARS also supports MARS melibea files to directly compare DL2 and DL3 data. It is also planned to add an interface to `pyIRF <https://github.com/cta-observatory/pyirf>`_ and `gammaPy <https://gammapy.org/>`_ at least for the CTLearn chain (MARS anaylsis should be converted via the DL3 converter). 

----

* Code : https://github.com/gae-ucm/pymars
* Author contact: Tjark Miener - tmiener@ucm.es
* License: GPL-3.0

----

Clone Repository with Git
^^^^^^^^^^^^^^^^^^^^^^^^^

Clone the pyMARS repository:

.. code-block:: bash

   cd </installation/path>
   git clone https://github.com/gae-ucm/pymars

Install Package with Anaconda/pypi
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Next, download and install `Anaconda <https://www.anaconda.com/download/>`_\ , or, for a minimal installation, `Miniconda <https://conda.io/miniconda.html>`_. Create a new conda environment that includes all the dependencies for LikelihoodCombiner:

.. code-block:: bash

   conda env create -f </installation/path>/pymars/environment.yml

Finally, install pyMARS into the new conda environment with pip:

.. code-block:: bash

   conda activate pymars
   cd </installation/path>/pymars
   pip install --upgrade .

(Under construction) Or install pyMARS via pypi (tested for Linux users):

.. code-block:: bash

   conda activate pymars
   pip install pymars

NOTE for developers: If you wish to fork/clone the respository and make changes to any of the pyMARS modules, the package must be reinstalled for the changes to take effect.

Installing as a conda package
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

(Under construction) To install it as a conda package, first install Anaconda by following the instructions here: https://www.anaconda.com/distribution/.

Then, create and enter a new Python3 environment with:

.. code-block:: bash

   conda create -n [ENVIRONMENT_NAME]
   conda activate [ENVIRONMENT_NAME]

From the environment, add the necessary channels for all dependencies:

.. code-block:: bash

   conda config --add channels conda-forge
   conda config --add channels menpo

Install the package:

.. code-block:: bash

   conda install -c tmiener pymars

This should automatically install all dependencies (NOTE: this may take some time, as by default MKL is included as a dependency of NumPy and it is very large).


Dependencies
^^^^^^^^^^^^

* astroPy
* ctaplot
* pyIRF (soon)
* gammaPy (soon)
* Python3
* Jupyter
* NumPy
* SciPy
* Pandas
* PyTables
* Matplotlib

Run pyMARS
----------

Run pyMARS from the command line:

.. code-block:: bash

   pymars-podie -h 


Uninstall pyMARS
----------------

Remove Anaconda Environment
^^^^^^^^^^^^^^^^^^^^^^^^^^^

First, remove the conda environment in which LikelihoodCombiner is installed and all its dependencies:

.. code-block:: bash

   conda remove --name pymars --all

Remove pyMARS
^^^^^^^^^^^^^

Next, completely remove pyMARS from your system:

.. code-block:: bash

   rm -rf </installation/path>/pymars
