.. figure::  _static/pyocto_logo_outlined.svg
   :align:   center

|

Welcome to PyOcto!
==================================

PyOcto is a high-throughput seismic phase associator.
The best way to get started with PyOcto is through our interactive examples.

.. raw:: html

    <embed>
        PyOcto basics:
        <a href="https://colab.research.google.com/github/yetinam/pyocto/blob/main/examples/01_basics.ipynb">
            <img src="https://colab.research.google.com/assets/colab-badge.svg" alt="Open In Colab"/>
        </a>
    </embed>
    <br>


.. raw:: html

    <embed>
        Velocity models:
        <a href="https://colab.research.google.com/github/yetinam/pyocto/blob/main/examples/02_velocity_models.ipynb">
            <img src="https://colab.research.google.com/assets/colab-badge.svg" alt="Open In Colab"/>
        </a>
    </embed>
    <br>


.. raw:: html

    <embed>
        Interfaces:
        <a href="https://colab.research.google.com/github/yetinam/pyocto/blob/main/examples/03_interfaces.ipynb">
            <img src="https://colab.research.google.com/assets/colab-badge.svg" alt="Open In Colab"/>
        </a>
    </embed>
    <br>

If you're trying to figure out how to set the configuration parameters,
check out our :ref:`guide on parameter choices<parameters>`.
For further details, please refer to the comprehensive documentation of
:ref:`the associator<associator>`, :ref:`the velocity models<velocity_models>`,
and :ref:`the data formats<data_formats>`.

You can install PyOcto directly through pip: ::

    pip install pyocto

Alternatively, you can install PyOcto from source by cloning the
`PyOcto Github repository <https://github.com/yetinam/pyocto>`_.

.. warning::

    PyOcto uses POSIX threads for threading. As these are not available on Windows,
    the Windows version is single-threaded. Therefore, we do not recommend running
    larger computations on Windows.

.. admonition:: Citation

   If you're using PyOcto in your work, please cite

      Münchmeyer, J. (2023). PyOcto: A high-throughput seismic phase associator. Preprint at `arxiv.org/abs/2310.11157 <https://arxiv.org/abs/2310.11157>`_.

.. toctree::
    :maxdepth: 2
    :caption: Contents:
    :hidden:

    self
    pages/associator.rst
    pages/velocity_models.rst
    pages/data_formats.rst
    pages/parameters.rst


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
