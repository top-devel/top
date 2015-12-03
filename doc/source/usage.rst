.. _usage:

#####
Usage
#####

TOP consists of a language that defines equations to be solved (see
:ref:`equation files<equation>`), a compiler wrapper designed to compile these
equation files, and a python module to actually run computations.

Compiling an Equation File
==========================

Compiling an equation file is performed by the compiler wrapper ``top-build``.
Equation files are attached to a given star model, you have to provide
``top-build`` with the model corresponding to the equation file with the
``--model=`` option.

**Example: poly_ester**

.. code-block:: shell

   top-build --model=poly_ester eq_poly_ester


Running TOP
===========

In order to use your newly compiled equation file, you need to use the top python
module:

.. code-block:: python

   import top                      # imports the top module 
   import matplotlib.pyplot as plt # import matplotlib for plots

   top.load('eq_poly_ester')       # loads your compiled equation file
   top.read_dati('dati')           # reads the parameter file `dati`
   model = top.dati.dirmodel       # dirmodel has to be an input in eq_file
   shift = top.dati.shift

   top.init_model(model)           # initializes the model stored in the directory
   top.run_arncheb(shift)          # runs Arnoldi-Chebyshev method

   # get the solutions in the first domain:
   for sol in range(0, top.get_nsol()): # for all solutions
       for var in top.get_vars(1):      # for each variable (in the fist domain)
           # get the solution number sol for variable named var in the first domain
           valp, vecp = top.get_sol(1, sol, var)
           plt.plot(vecp)
           plt.title(var)
           plt.show()
           print('Variable: %s, valp: %e' % (var, valp))
