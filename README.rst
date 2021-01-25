.. code::

       ______      _____ ____ 
      / ____/___ _/ ___// __ \
     / /   / __ `/\__ \/ / / /
    / /___/ /_/ /___/ / /_/ / 
    \____/\__,_//____/\___\_\ 

|pipeline status| |coverage report| |black| |rtd| |gpl|

|pypi-version| |pypi-python| |pypi-wheel| |pypi-downloads| |deps|

.. |pipeline status| image:: https://gitlab.inria.fr/soliman/casq/badges/master/pipeline.svg
   :target: https://gitlab.inria.fr/soliman/casq/commits/master
   :alt: pipeline status

.. |coverage report| image:: https://gitlab.inria.fr/soliman/casq/badges/master/coverage.svg
   :target: https://gitlab.inria.fr/soliman/casq/commits/master
   :alt: coverage report

.. |black| image:: https://img.shields.io/badge/code%20style-black-000000.svg
   :target: https://github.com/python/black
   :alt: Code style: black

.. |rtd| image:: https://readthedocs.org/projects/casq/badge/?version=latest
   :target: https://casq.readthedocs.io/en/latest/?badge=latest
   :alt: Documentation Status

.. |gpl| image:: https://img.shields.io/pypi/l/casq
   :target: https://gitlab.inria.fr/soliman/casq/raw/master/LICENSE
   :alt: PyPI - License

.. |pypi-version| image:: https://img.shields.io/pypi/v/casq
   :target: https://pypi.org/project/casq/
   :alt: PyPI

.. |pypi-python| image:: https://img.shields.io/pypi/pyversions/casq
   :alt: PyPI - Python Version
   :target: https://pypi.org/project/casq/

.. |pypi-wheel| image:: https://img.shields.io/pypi/wheel/casq
   :target: https://pypi.org/project/casq/
   :alt: PyPI - Wheel

.. |pypi-downloads| image:: https://img.shields.io/pypi/dm/casq
   :target: https://pypi.org/project/casq/
   :alt: PyPI - Downloads

.. |deps| image:: https://img.shields.io/librariesio/release/pypi/casq
   :target: https://pypi.org/project/casq/
   :alt: Libraries.io dependency status for latest release

**CaSQ** converts `CellDesigner`_ models to Boolean models encoded in
`SBML-Qual`_ with a rather strict semantics defined in a
*submitted* article.

.. _`CellDesigner`: http://celldesigner.org
.. _`SBML-Qual`: http://sbml.org

Install
=======

CaSQ is provided as a Python (>= 3.5) package, you can install it from the `Python package index`_ with ``pip``, ``conda`` or your Python package manager of choice:

.. _`Python package index`: https://pypi.org/project/casq/

.. code:: bash

   $ pip3 install casq

Usage
=====

Just follow the instructions::

   $ casq --help
   usage: casq [-h] [-v] [-d] [-c] [-s] [-r S] [infile] [outfile]

   Convert CellDesigner models to SBML-qual with a rather strict semantics.
   Copyright (C) 2019-2021 Sylvain.Soliman@inria.fr GPLv3

   positional arguments:
     infile            CellDesigner File
     outfile           SBML-Qual File

   optional arguments:
     -h, --help        show this help message and exit
     -v, --version     show program's version number and exit
     -d, --debug       Display a lot of debug information
     -c, --csv         Store the species information in a separate CSV file
     -s, --sif         Store the influence information in a separate SIF file
     -r S, --remove S  Delete connected components in the resulting model if
                       their size is smaller than S. A negative S leads to keep
                       only the biggest(s) connected component(s)
