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

   import top                       # imports the top module 
   import numpy as np               # imports numpy

   p = top.load('eq_poly_ester')    # loads your compiled equation file
   p.read_dati('dati')              # reads the parameter file `dati`

   model = 'model/'                 # path to the model
   shift = p.dati.shift

   m = p.init_model(model)          # initializes the model stored in the directory
   r = p.run_arncheb(shift)         # runs Arnoldi-Chebyshev method

   # get the solutions
   for i in range(0, r.nsol):       # for all solutions
       for var in p.get_vars(0):    # for each variable (in the fist domain)
           # plot the solution
           r.plot(0, i, v)          # quick plot of the solution

           # plot an expression of the solution:

           # get the solution
           val, vec, l = r.get_sol(0, i, v)

           # get the grid
           radius, theta = r.get_grid()
           cost = np.cos(theta)

           # get a field from the model
           h = m['hh']

           # project the solution onto the grid
           gv = top.leg.eval2d(vec, cost, l[0], 2, r.dati['m'])

           # actually plot the expression
           r.plot_val(gv*np.sqrt(h)**r.dati['pindex'])



For a detailed description of all functionalities available through the python
module, see :ref:`python API<api>`.
