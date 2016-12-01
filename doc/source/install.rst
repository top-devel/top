#######
Install
#######

Quick Install Guide
===================

TOP uses the standard ``autotools`` (autoconf_, automake_) install procedure,
you should be able to install it using:

.. _autoconf: https://www.gnu.org/software/autoconf/
.. _automake: https://www.gnu.org/software/automake/

.. code-block:: shell

    # get the code:
    git clone https://gitlab.com/top-dev/top.git

    # enter the source code directory:
    cd top

    # download the optional dependencies:
    git submodule init
    git submodule update

    # prepare the configure script:
    ./bootstrap

    # configure and install TOP:
    ./configure && make install

You can get more details and options reading the following sections:
:ref:`getting the code<download>`, :ref:`configure<configure>` and
:ref:`compile and install<build>`.


.. _download:

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

This will download the latest version of the code in a directory named ``top``.
Enter this directory:

.. code-block:: shell

    cd top

You can optionally download the parser to write oscillation equations with the
new equation format:

.. code-block:: shell

    git submodule init
    git submodule update

Run the ``bootstrap`` script that will create the configure script:

.. code-block:: shell

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
In order to install TOP, you will need:

* a Fortran compiler supporting procedure interface (``gfortran (>=4.9)``)
* a recent version of Python (``python (>=2.7)``)
* the program ``f2py``, usually shipped with ``numpy``
* the following python modules: ``numpy`` and ``h5py``

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
* ``PYTHON``: the python interpreter to use (e.g. ``PYTHON=python3``)
* ``--prefix=``: this option allows you to set TOP's install directory (by
  default the prefix is set to ``$HOME/local``)


**Example:**

If you want to use Intel compiler (:samp:`ifort`) and the ATLAS library
(installed in ``$HOME/local/lib``), you want to configure with the following
command line:

.. code-block:: shell

   ./configure FC=ifort LDFLAGS=-L$HOME/local/lib LIBS=-ltatlas

.. _build:

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
See :ref:`usage<usage>`.

Using ``libester``
==================

In order to use ESTER_ stellar models, TOP needs to find where ESTER was
installed on your system.
In order to tell TOP's configure script where to find ``libester``, you need to
provide it with the options:
``LDFLAGS=-L$PATH_TO_ESTER/lib`` and ``CPPFLAGS=-I$PATH_TO_ESTER/include``.

For instance if ESTER was installed in ``$HOME/local``, TOP should be able to
find it if you configure with:

.. code-block:: shell

    ./configure LDFLAGS=-L$HOME/local/lib CPPFLAGS=-I$HOME/local/include

.. _ESTER: http://ester-project.github.io/ester/
