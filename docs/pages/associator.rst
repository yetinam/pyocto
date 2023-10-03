.. _associator:

The OctoAssociator class
===========================

A PyOcto workflow consists of two steps. First, you need to create an associator instance.
In this step, you set the configuration for the model.
Second, you run the associator by providing it with a list of picks and station metadata.
If you've previously used the GaMMA or REAL associators, you might be interested in
using the compatibility interfaces for this step.

.. autoclass:: pyocto.associator.OctoAssociator
    :members:
    :undoc-members:
    :show-inheritance:
