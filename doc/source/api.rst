.. _api:

Python API
==========

.. autoclass:: top.load

   .. automethod:: top.load.get_version
   .. automethod:: top.load.read_dati
   .. automethod:: top.load.init_model
   .. automethod:: top.load.run_arncheb
   .. automethod:: top.load.get_grid
   .. automethod:: top.load.get_zeta
   .. automethod:: top.load.get_sol
   .. automethod:: top.load.get_vars
   .. automethod:: top.load.get_results
   .. automethod:: top.load.write_output


.. autoclass:: top.results

   .. automethod:: top.results.get_sol
   .. automethod:: top.results.get_grid
   .. automethod:: top.results.save
   .. automethod:: top.results.plot_val
   .. automethod:: top.results.plot
   .. automethod:: top.results.read_model


.. _model:

Star Models
-----------

In order to access fields of the model, you can use square bracket operator `[]`
with the name of the field between brackets (e.g, `m['h']` would return the
enthalpy of the model `m`).

But you need to know what are the fields defined for each model. See
:ref:`models<models>` for a list of supported models in TOP, and a list a fields defined by each model.

.. autoclass:: top.model
