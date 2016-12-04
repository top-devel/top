.. _api:

Python API
==========

Equation
--------

The :class:`top.load` class in responsible for loading a previously compiled
equation file, loading or setting parameters, loading stellar models and solving
the underlying eigenvalue problem.

.. autoclass:: top.load
    :members:


Results Objects
---------------

.. autoclass:: top.results
    :members:

.. _model:

Star Models
-----------

In order to access fields of the model, you can use square bracket operator `[]`
with the name of the field between brackets (e.g, `m['h']` would return the
enthalpy of the model `m`).

But you need to know what are the fields defined for each model. See
:ref:`models<models>` for a list of supported models in TOP, and a list a fields defined by each model.

.. autoclass:: top.model
    :members:
