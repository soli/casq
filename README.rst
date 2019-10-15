.. code::

       ______      _____ ____ 
      / ____/___ _/ ___// __ \
     / /   / __ `/\__ \/ / / /
    / /___/ /_/ /___/ / /_/ / 
    \____/\__,_//____/\___\_\ 

|pipeline status| |coverage report| |black|

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

**CaSQ** converts `CellDesigner`_ models to Boolean models encoded in
`SBML-Qual`_ with a rather strict semantics defined in a
currently-being-written article.

.. _`CellDesigner`: http://celldesigner.org
.. _`SBML-Qual`: http://sbml.org

Install
=======

CaSQ is provided as a Python (>= 3.5) package, you can install it with ``pip``,
``conda`` or your Python package manager of choice:

.. code:: bash

   $ pip3 install casq

Usage
=====

Just follow the instructions:

.. code:: bash

   $ casq --help
