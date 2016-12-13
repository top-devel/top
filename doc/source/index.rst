###################
TOP's documentation
###################

Overview
########

TOP consists of a python module, a compiler and several *templates* (see
:ref:`models`) to read stellar models.
This allows users to write their own set of equations in a flexible way.

The basic workflow with TOP is the following:

1. write an equation file (see :ref:`equation`)
2. compile this file with ``top-build``
3. compute oscillations modes and frequencies with ``top`` python module (see
   :ref:`api`):

   - read input parameters
   - read a stellar model
   - run the Arnoldi-Chebyshev method

.. figure:: overview.png
    :width: 400px
    :align: center
    :figwidth: 500px

    TOP's Software Architecture


Documentation
#############

.. toctree::
   :maxdepth: 2

   install
   usage
   model
   equation
   api
   examples
   download

..   hacking

Credits
#######
This public version of the TOP code is part of the SPACEInn_ Project (Exploitation of Space Data for Innovative Helio- and Asteroseismology), initiated by the European Helio- and Asteroseismology Network (HELAS_) and financed by the European Union under the Seventh Framework Programme (FP7 project n° 312844). This work also benefits of the French Spatial Agency (CNES_) through its support to the *Plato* mission.

.. image:: spaceinn_logo.jpg
   :alt: SPACEInn
   :height: 100px
   :target: http://www.spaceinn.eu/
   :align: left
.. image:: logo-irap.png
   :alt: IRAP
   :height: 100px
   :target: http://www.irap.omp.eu/
   :align: left
.. image:: Logo-UPS.jpg
   :alt: UPS Toulouse
   :height: 100px
   :target: http://www.univ-tlse3.fr/
   :align: left
.. image:: logo_cnrs.png
   :alt: CNRS
   :height: 100px
   :target: http://www.cnrs.fr/
   :align: left

.. _SPACEInn: http://www.spaceinn.eu/
.. _HELAS: https://www.helas-project.eu/
.. _CNES: https://cnes.fr/en
