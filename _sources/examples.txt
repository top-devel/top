##############
Usage Examples
##############

Surface
#######

Description
===========

The directory example contains a subdirectory name ``surface`` that proposes python softwares and scripts:

* ``patch_model.py`` is a module patched an averaged 3-D atmosphere on a 1-D stellar structure

* ``run_model.py`` is an example of use of the ``patch_model`` module followed by the frequency computation for the unpatched and patched model. The example is based on the solar cases. The reference model is the model S by by Christensen-Dalsgaard, J. et al. (1996, Science 272, 1286), the stellar atmosphere has been computed by Belkacem, K.,Samadi , R.,Kupka, F., Grimm-Strele, H. (in preparation) with the ANTARES code Muthsam et al. (2010, New Ast. 15, 460). At the end of the example the script compared computed frequencies with GOLF observation (Courtesy R.A. Garcia, CEA-Saclay)

* ``script.sh`` is a script compiling top and running ``run_model.py``.

Results
=======

Here are some results obtained with ``run_model.py``
 
.. figure:: figure_0:rho.png
   :width: 500px
   :align: center
   :figwidth: 600px
   :alt: density

   Density profiles near the surface.

.. figure:: figure_10:patched_result.png
   :width: 500px
   :align: center
   :figwidth: 600px
   :alt: comparison

   The frequencies of radial modes observed with GOLF compared to the ones computed with patched and unpatched models.

These results are based on an approach similar to the one used by  Rosenthal et al. (1999, A&A 351, 689) or Piau et al. (2014, MNRAS 437, 164).

