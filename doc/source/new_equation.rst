.. _new_equation:

###################
New Equation Format
###################

Introduction
############

This equation format aims at simplifying the way to write oscillation equations.
With this format equation are written in there mathematical expression after
projection on the spherical harmonics.

Here is an example of such equation if the 1D case:

.. code::

    lh*(lh+1) * Et = 
        avg(r/Gamma1) * dP_P    +
        avg(r)        * Er'     +
        2             * Er
    with (r=1)
        dr(Er, -1) = lh * dr(Et, -1)
        at r = 0


.. note::

    TOP solves eigenvalue problems. Therefore equations written for TOP must be
    linear.

Language Description
####################

.. _param_def:

The Parameters Section
======================

At the beginning of an equation file, one is able to define several parameter
that can be used within the definition of equation. These parameters can be
defined with the ``input`` keyword, followed by its type (``double``, ``int``
or ``string``) ant its name:

:Example:

.. code::

    input double mass
    input double rota

Variable and Model Definition
=============================

In order to *understand* the meaning of the equation to appear, TOP's compiler
need to know what are the variables of the problem, and what are the field of the
model.

.. _var_def:

Variable Definition
-------------------

In order to define variables of the problem, one has to define them using the
``var`` keyword followed by a comma separated list of names.

.. note::

    2 or 3 variables surrounded by parenthesis means that they are component of
    vector.

:Example:

.. code::

    var Phi, PhiP, (Er, Et), dP_P

Will define 5 variables, ``Er`` and ``Et`` being first and second component
of a vector.


.. _model_def:

Model's Field definition
------------------------

In the same way variables can be defined, fields and scalar defined from the
stellar model can be defined with ``field`` or ``scalar`` keywords:

:Example:

.. code::

    field pm, g_m, r, rhom, dg_m, rhom_z
    scalar Gamma1, Lambda

Equation Definition
===================

After the definition/declaration section, we can start defining the equation
after the ``in`` keyword.

Defining an equation
--------------------

In order to add an equation in the system, one can use the ``equation`` keyword,
followed by the name of the equation, followed by a ``:`` (colon) and the
expression of the equation.

:Example:

.. code::

    equation eqdP_P:
    lh*(lh+1) * Et = 
        avg(r/Gamma1) * dP_P    +
        avg(r)        * Er'     +
        2             * Er

.. note::

    Every identifier (*i.e.,* name) involved in an equation need to be defined
    (either as a :ref:`variable<var_def>`, a :ref:`field or a scalar<model_def>` or
    :ref:`a parameter<param_def>`).

Boundary Condition
------------------

In order to define boundary condition, one simply need to define of after the
equation with the following syntax:

:Syntax:
    ``with (r=numerical_location) epxression at r = location`` where:

    :numerical_location: is the line of the matrix to be replaced with the
             boundary condition.
    :expression: if the expression of the boundary condition.
    :location: is the *physical* location of the boundary condition.

:Example:

.. code::

    equation eqdP_P:
    lh*(lh+1) * Et = 
        avg(r/Gamma1) * dP_P    +
        avg(r)        * Er'     +
        2             * Er
    with (r=1)
        dr(Er, -1) = lh * dr(Et, -1)
        at r = 0

Internal variables, and Functions
=================================

A few functions and variables are already defined with TOP and can be used
without prior declarations, here is a list of such symbols:

:``dr``:

    :Syntax: ``dr(var, order)``

    :Semantics: derivative of ``var`` of order ``order``:

    :Example:
        ``dr(Phi, 2)``


        stands for :math:`\frac{\partial^2\Phi}{\partial r^2}`


    .. note::

        Radial derivatives can also be expressed with the ``'`` (apostrophe)
        post-fixed operator: ``dr(Phi, 2)`` and ``Phi''`` are two notations strictly
        equivalent.

--------------------------------------------------------------------------------

:``fp``:

    :Syntax: ``fp``

    :Semantics: ``fp`` is the eigenvalue of the problem. It should appear in
        equation definition

    :Example: ``fp^2 * r * Et - pm/rhom * dP_P - Phi - g_m = 0``

--------------------------------------------------------------------------------

:``avg``:

    :Syntax: ``avg(expr)``

    :Semantics: average of expression ``expr`` on the point of the grid used for for
                derivation or interpolation (therefore it depends on the numerical
                scheme used).

    :Example: ``avg(r/Gamma1) * dP_P``

Comments
========

Comments can be added in equation file using a pound sign (``#``), the remaining
of the line will be ignored.

:Example:

    .. code::

        # define the first equation
        equation eqdP_P:
        lh*(lh+1) * Et =            # this is the LHS of the equation
            avg(r/Gamma1) * dP_P +  # this the RHS
            avg(r) * Er'
