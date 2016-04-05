.. _models:

##################
Star Models in TOP
##################

This page describes the star models supported by TOP as well as the fields
defined in the model.

Polytropic Model
================

**Name:** ``poly_ester``

Fields:
 * ``hh``
 * ``hht``
 * ``hhz``
 * ``hhzz``
 * ``hhzt``
 * ``lnhht``
 * ``lambda``
 * ``alpha``
 * ``aplat``
 * ``omega_K``
 * ``omga``
 * ``r_t``
 * ``r_z``
 * ``r_map``
 * ``re_t``
 * ``re_z``
 * ``re_map``
 * ``r_zz``
 * ``r_zt``
 * ``r_tt``
 * ``re_zz``
 * ``re_zt``
 * ``re_tt``
 * ``zeta``
 * ``cost``
 * ``sint``
 * ``cott``

ESTER Model
===========

**Name**: ``ester``

Fields:
 * ``zeta``
 * ``theta``

CESAM Model
===========

**Name**: ``cesam``

Fields:
 * ``rhom``: density :math:`\rho`
 * ``rhom_z``: :math:`\frac{\partial{\rho}}{\partial{\zeta}}`
 * ``rhom_t``: :math:`\frac{\partial{\rho}}{\partial{\theta}}`
 * ``pm``: pressure :math:`p`
 * ``pm_z``: :math:`\frac{\partial{p}}{\partial{\zeta}}`
 * ``pm_t``: :math:`\frac{\partial{p}}{\partial{\theta}}`
 * ``Gamma1``: first adiabatic exponent
 * ``NN``: Brunt-Vaisala frequency
 * ``NNr``: :math:`\frac{\partial{NN}}{\partial{\zeta}}`
 * ``NNt``: :math:`\frac{\partial{NN}}{\partial{\theta}}`
 * ``pe``: pressure potential
 * ``pe_z``: :math:`\frac{\partial{pe}}{\partial{\zeta}}`
 * ``pe_t``: :math:`\frac{\partial{pe}}{\partial{\theta}}`
