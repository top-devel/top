.. _equation:

##############
Equation Files
##############

There are 2 different formats to write oscillation equations in TOP: the *legacy*
one and the *newer* one.

Legacy Equation Format
======================

In this format, equations consists of series of commands that build the
numerical system.
In the following, we describe theses commands:

input
-----
This command allow the user to define a parameter for the set of equation to be
written.

:Syntax:
    ``input parameter format``
:Arguments:
    - ``parameter``: name of the parameter.
    - ``format``: FORTRAN format specifier. This is used to read the
      parameter input file.
:Example:
    ``input mass 0pf5.2``

stamp
-----
Used to define a string to appear in the output files (deprecated).

:Syntax:
    ``stamp string``
:Arguments:
    - ``string``: string to appear in output files.
:Example:
    ``stamp eq_ESTER_all_lagrange``

definition
----------
Allow users to define constants

:Syntax:
    ``definition type name value``
:Arguments:
    - ``type``: type of the constant (``integer``, ``double_precision``,
      ``complex``)
    - ``name``: name of the constant
    - ``value``: value of the constant
:Example:
    ``definition double_precision gamma_p 1d0 + 1d0/pindex``
    will define ``gamma_p`` with a value of :math:`1+\frac{1}{pindex}`
    (where :math:`pindex` has to be another variable (``definition`` or
    ``input``) 

eqlist
------
Defines the name of equation to be defined in the equation file

:Syntax:
    ``eqlist eq1 eq2 eq3 ... # and so on``
:Arguments:
    - ``eq1``: name of the equation
:Example:
    ``eqlist eqEr eqdP eqPhi eqPhiP``
    defines 4 equations named eqEr, eqdP, eqPhi and epPhiP

varlist
-------
Defines variables of the equation set

:Syntax:
    ``varlist var1 var2 var3 ... # and so on``
:Arguments:
    - ``var1``: name of the variable
:Example:
    ``varlist Er dP Phi PhiP``
    defines 4 variables named Er, dP, Phi and PhiP

leq
-------
In TOP equation are projected into the spherical harmonic basis. This command is
use to define the starting :math:`l` for this projection.

:Syntax:
    ``leq eqName value``
:Arguments:
    - ``eqName``: name of the equation
    - ``value``: starting value of :math:`l`
:Example:
    ``leq eqEr abs(m)+iparity``

lvar
-------
In TOP variables are projected into the spherical harmonic basis. This command is
use to define the starting :math:`l` for this projection.

:Syntax:
    ``lvar varName value``
:Arguments:
    - ``varName``: name of the variable
    - ``value``: starting value of :math:`l`
:Example:
    ``lvar Er abs(m)+iparity``

equation
--------
This command is used to start defining an equation. This means that further
command in the equation file will apply to the *current* equation.

:Syntax:
    ``equation eqName``
:Arguments:
    - ``eqName``: name of the equation
:Example:
    ``equation eqEr``

sub
--------
This is use to insert a term in the *current* equation: this term will be
computed by calling a FORTRAN subroutine.

:Syntax:
    ``sub type power routine variable``
:Arguments:
    - ``type``: the type of term see :ref:`type of terms in TOP<terms>`.
    - ``power``: the power of the eigenvalue preceded by a ``w``.
    - ``routine``: name of the FORTRAN subroutine to be called to compute the
      coupling coefficient.
    - ``variable``: name of the variable involved in the coupling. Further
      characters can be used indicate radial derives. For instance,
      ``Er'`` mean :math:`\frac{\partial Er}{\partial r}`. Higher derivative
      order can be achieved either by chaining the ``'`` character or with the
      ``^`` character followed by the derivative order: ``Er^2`` is equivalent
      to ``Er''``.
:Example:
    ``sub rtt w1 Illm(sint/roz, $a, $leq, $lvar) u``: this basically add the
    term :math:`\omega \iint(\frac{sin(\theta)}{roz}) * u` in the current
    equation, where :math:`\omega` is the eigenvalue.


--------------------------------------------------------------------------------

.. _terms:

:Type of terms in TOP:
    - ``s``: scalar term
    - ``r``: term only depending on :math:`r` (or :math:`\zeta`) (only available for
      1D equations)
    - ``rt``: term depending on :math:`r` (or :math:`\zeta`) and :math:`l` (only in
      2D equations)
    - ``rtt``: term depending on :math:`r` (or :math:`\zeta`), :math:`l` and
      :math:`l'` (only available in 2D equations)

New Equation Format
===================
