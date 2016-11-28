######################
Legacy Equation Format
######################

In this format, equations consists of series of commands that build the
numerical system.

Overview
========

In order to write oscillation equations with this format, one need to write them
projected onto the spherical harmonics basis.

:Equation example:

.. math::

    \lambda b_m^l = \sum_{l'=|m|}^{\infty}
    - \iint_{4\pi} \{Y_l^m\}^*Y_{l'}^m d\Omega \partial_\zeta u_m^{l'}
    - \iint_{4\pi} \frac{2 \zeta H + \zeta^2 N H_\zeta}{r^2 r_\zeta} \{Y_l^m\}^*Y_{l'}^m d\Omega u_m^{l'}
    + ...

In order to write such an equation in TOP, one need to split it into terms and
add them incrementally in the numerical system. The purpose of TOP, in
particular ``top-build`` is to ease writing such equations.
Several commands can help inserting such terms in the numerical system, see
section :ref:`commands` for a list of commands available in TOP's language.

Few other features can help understanding how to write equations in TOP: for a
description of the type of terms TOP can handle, see section :ref:`terms`. And
see section :ref:`preprocessing` for a description of string pre-processing and
special variables.

.. _terms:

Type of Terms
-------------

For each term to be added in the numerical system, one need to provide TOP with
its type. The type of a term specify whether it depends on nothing (scalar), the
radial coordinate :math:`r`, :math:`l` or :math:`l'`.

The type of terms currently supported by TOP are:

:``s``: scalar term
:``r``: term only depending on :math:`r` (or :math:`\zeta`) (only available for
        1D equations)
:``tt``: term depending on :math:`l` and :math:`l'` (only available in 2D
         equations)
:``rt``: term depending on :math:`r` (or :math:`\zeta`) and :math:`l` (only in
         2D equations)
:``rtt``: term depending on :math:`r` (or :math:`\zeta`), :math:`l` and
          :math:`l'` (only available in 2D equations)

.. _preprocessing:

String pre-processing 
=====================

:``$a``:
    stands for the current coupling matrix. This variable can only be used in a
    term definition

:``$prev``:
    the previous coupling matrix. It correspond to variable ``$a`` of the last
    term definition.

:``$leq``:
    the array containing :math:`l` values of the current equation.

:``$lvar``:
    the array containing :math:`l` values of the current variable.

:``$leq(x)``:
    is the :math:`x^{th}` :math:`l` value of the current equation.

:``$lvar(x)``:
    is the :math:`x^{th}` :math:`l` value of the current variable.

:``$eq``:
    index of the current equation.

:``$var``:
    index of the current variable.

:``$nr``:
    radial resolution of the problem.

:``$i``:
    indices of the radial coordinate (this will result in the generation of a
    FORTRAN loop over all radial points).

:``$j1``:
    indices of the horizontal coordinate :math:`l` (this will result in the
    generation of a FORTRAN loop over all values of :math:`l`).

:``$j2``:
    indices of the horizontal coordinate :math:`l'` (this will result in the
    generation of a FORTRAN loop over all values of :math:`l'`).

.. _commands:

Commands
========

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
Used to define a string to appear in the output files.

:Syntax:
    ``stamp string``
:Arguments:
    - ``string``: string to appear in output files.
:Example:
    ``stamp eq_ESTER_all_lagrange``

definition
----------
Used to define named constants.

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
Defines the name of equation to be defined in the equation file.

:Syntax:
    ``eqlist eq1 eq2 eq3 ... # and so on``
:Arguments:
    - ``eq1``: name of the equation
:Example:
    ``eqlist eqEr eqdP eqPhi eqPhiP``
    defines 4 equations named eqEr, eqdP, eqPhi and epPhiP

varlist
-------
Defines variables of the equation set.

:Syntax:
    ``varlist var1 var2 var3 ... # and so on``
:Arguments:
    - ``var1``: name of the variable
:Example:
    ``varlist Er dP Phi PhiP``
    defines 4 variables named Er, dP, Phi and PhiP

leq
...
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
----
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
...
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

subbc
-----
This is use to insert a boundary condition term in the *current* equation: this
term will be computed by calling a FORTRAN subroutine.

:Syntax:
    ``subbc type location power routine variable(index)``
:Arguments:
    - ``type``: the type of term see :ref:`type of terms in TOP<terms>`.
    - ``location``: the location where the boundary condition should be inserted
      in the **numerical** system. This is basically tells the line in the
      matrix to be replaced with the boundary condition.
    - ``power``: the power of the eigenvalue preceded by a ``w``.
    - ``routine``: name of the FORTRAN subroutine to be called to compute the
      coupling coefficient.
    - ``variable``: name of the variable involved in the coupling. Further
      characters can be used indicate radial derives. For instance,
      ``Er'`` mean :math:`\frac{\partial Er}{\partial r}`. Higher derivative
      order can be achieved either by chaining the ``'`` character or with the
      ``^`` character followed by the derivative order: ``Er^2`` is equivalent
      to ``Er''``.
    - ``index'`` radial coordinate of the boundary condition.
:Example:
    ``subbc tt nr w0 Illmbc(hhz(1, :), $a, $leq, $lvar) v(1)``, here we can see
    that ``location`` and ``index`` are different: the boundary condition is
    imposed at the center (``v(1)`` stands for :math:`v` at :math:`r=0`), but in
    the **numerical** system, the condition is imposed on the last line of the
    matrix.

term
----
Used to insert a term in the equation.

:Syntax:
    ``term type power expression  variable``
:Arguments:
    - ``type``: the type of term see :ref:`type of terms in TOP<terms>`.
    - ``power``: the power of the eigenvalue preceded by a ``w``.
    - ``expression``: the mathematical expression of the term to be inserted.
    - ``variable``: name of the variable involved in the coupling. Further
      characters can be used indicate radial derives. For instance,
      ``Er'`` mean :math:`\frac{\partial Er}{\partial r}`. Higher derivative
      order can be achieved either by chaining the ``'`` character or with the
      ``^`` character followed by the derivative order: ``Er^2`` is equivalent
      to ``Er''``.
:Example:
    ``term s w0 -2d0 Pi''``: this would insert the term
    :math:`-2\frac{\partial^2 Pi}{\partial r^2}` in the current equation.

termbc
------
Used to insert a term in a boundary condition of the system.

:Syntax:
    ``termbc type location power expression  variable(index)``
:Arguments:
    - ``type``: the type of term see :ref:`type of terms in TOP<terms>`.
    - ``power``: the power of the eigenvalue preceded by a ``w``.
    - ``location``: the location where the boundary condition should be inserted
      in the **numerical** system. This is basically tells the line in the
      matrix to be replaced with the boundary condition.
    - ``expression``: the mathematical expression of the term to be inserted.
    - ``variable``: name of the variable involved in the coupling. Further
      characters can be used indicate radial derives. For instance,
      ``Er'`` mean :math:`\frac{\partial Er}{\partial r}`. Higher derivative
      order can be achieved either by chaining the ``'`` character or with the
      ``^`` character followed by the derivative order: ``Er^2`` is equivalent
      to ``Er''``.
    - ``index'`` radial coordinate of the boundary condition.
:Example:
    ``termbc t $nr w0 1d0 Phi($nr)'``: this would insert the term
    :math:`\Phi(r=surf)` in the boundary condition. (The last line of the
    matrix would be replaced with this boundary condition).

instruction
-----------
Used to add ad-hoc FORTRAN instruction in the module responsible for computing
coupling integrals.

:Syntax:
    ``instruction fortran``
:Arguments:
    - ``fortran``: the FORTRAN instruction to be inserted.
:Example:
    ``instruction call modify_l0($prev, $nr, abs(m)+iparity)``: will insert the
    code ``call modify_l0(dm(1)%artt(:, :, :), grd(1)%nr, abs(m)+iparity)``.
    See :ref:`preprocessing`.
