.. _getting_started-Installation:


============
Installation
============

The following sections describe how to install the scflow framework.

------------------
Conda Installation
------------------

The preferred method for installation is through conda. Currently this installation is still in working progress. Preferably the installation should be in a separate environment.::

  # The conda environment is currently not ready but we are working on this
  #conda create -n scflow -c cgat scflow
  #conda activate scflow
  #scflow --help

----------------
Pip installation
----------------

You can install scflow using pip, this will only install the package without any dependancies, which will have to be installed separately.::

  pip install scflow

-------------------
Manual Installation
-------------------

The repository can also be installed manually, but dependancies will need to be installed separately.::

  git clone https://github.com/Acribbs/scflow.git
  python setup.py install
  scflow --help

-----------------------------
Installing additional software
-----------------------------

We always recommend using conda to install additional software where possible.

This can be easily performed by::

  conda search <package>
  conda install <package>
