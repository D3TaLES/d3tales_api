
Welcome to the |D3TaLES| API documentation!
===========================================

.. |D3TaLES| replace:: D\ :sup:`3`\ TaLES
.. image:: D3TaLES_api_logo.png
   :target: D3TaLES_api_logo.png
   :alt: D3TaLES API logo


API designed for the `D3TaLES project <https://d3tales.as.uky.edu/>`_. Includes data processing, database access functions, and property calculators.

Modules
-------

This API contains three modules:


* ``D3database``\ : Interface for interacting with the |D3TaLES| databases. Includes
  queries, data insertion with validation, and a REST API interface.

* ``Processors``\ : Processing code for extracting data from computational and experimental
  data files as well as `ExpFlow <https://d3tales.as.uky.edu/expflow>`_ run objects.

* ``Calculators``\ : Property calculators and plotters for cyclic voltammetry


.. toctree::
   :maxdepth: 3
   :caption: Contents:

   installation
   D3database
   Processors
   Calculators
   d3tales_api


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

