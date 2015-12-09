#######
Install
#######

Getting the code
================

From the git repository
-----------------------
The best way to get the latest version of the code is to download it from the
repository:

.. code-block:: shell

   # with ssh:
   git clone git@gitlab.com:top-dev/top.git

   # with https:
   git clone https://gitlab.com/top-dev/top.git

This will download the last version of the code in a directory named ``top``.
Enter this directory and create the configure script by running the
``bootstrap``
script.

.. code-block:: shell

   cd top
   ./bootstrap

.. note::

   In order to run the ``bootstrap`` script you will need to have ``autoconf``
   (>=2.59), ``automake`` (>=1.9) and ``libtool`` installed.


You can then proceed to the :ref:`configure<configure>` steps.

From a source archive
---------------------
If you don't need the latest version of TOP, you can use a source archive.
Extract your source archive and enter to source directory:

.. code-block:: shell

   tar xvjz top-x.y.tar.bz2
   cd top-x.y

And proceed to the :ref:`configure<configure>` steps.


.. _configure:

Configure
=========

**Prerequisites:**

The configure script allows you to configure the build environment of TOP.
The only prerequisite to install TOP is to have a BLAS and a LAPACK library
installed.
The program ``f2py``, usually shipped with ``numpy``, is also required.

Configure will try to detect the libraries installed in your system, if it fails
to find both a BLAS and a LAPACK library it will return an error.
You can try to re-run configure with some of the following option to help it
find you libraries:

**Configure options:**

* ``FC``: allows you to choose your Fortran compiler (e.g. ``FC=gfortran``)
* ``LDFLAGS``: sets linker flags. This can be used to specify libraries search
  directory (e.g. ``LDFLAGS=-L$HOME/local/lib``)
* ``LIBS``: what libraries should be linked with TOP. (e.g. ``LIBS=-ltatlas``)
* ``CPPFLAGS``: preprocessor flags, this can be used to tell the compiler where
  to find header files (e.g. ``CPPFLAGS=-I$HOME/local/include``)
* ``--prefix=``: this option allows you to set TOP's install directory (by
  default the prefix is set to ``$HOME/local``)


**Example:**

If you want to use Intel compiler (:samp:`ifort`) and the ATLAS library
(installed in ``$HOME/local/lib``), you want to configure with the following
command line:

.. code-block:: shell

   ./configure FC=ifort LDFLAGS=-L$HOME/local/lib LIBS=-ltatlas

Compile & Install
=================

After running successfully the configure script, you can compile and install TOP by running:

.. code-block:: shell

   make install

TOP is composed of a compiler wrapper ``top-build`` installed in
``$prefix/bin``, a few libraries installed in ``$prefix/lib`` and a python
module installed in ``$prefix/lib/python-version/site-packages/top``.

As few examples are also availiable in ``$prefix/share/top/models``

.. note::

   You can source the shell script ``activate-top.sh`` created in the directory
   where you compiled TOP to set up the environment variables PATH,
   LD_LIBRARY_PATH and PYTHONPATH with the path where TOP was installed.

Check you Install
=================
Seed :ref:`usage<usage>`.

