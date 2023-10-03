.. _velocity_models:

Velocity models
=========================

PyOcto integrates two velocity models, homogeneous and 1D layered.
As a rule of thumb, the homogeneous model is absolutely sufficient for
regions with only shallow (<30 km) seismicity and of limited size (<250 km diameter).
Here the approximation of constant velocity is usually good enough.
When working in zones with deeper seismicity, e.g., subduction zones,
or with a larger seismic network, the 1D model usually performs better.
Nonetheless, the homogeneous model usually still is sufficient for a
decent catalog.

.. warning::
    PyOcto uses a local coordinate transform to calculate distances.
    For very large study areas (> 2000 km diameter), this might lead to biases in distance calculation.

.. autoclass:: pyocto.associator.VelocityModel0D
    :members:
    :undoc-members:
    :show-inheritance:

.. autoclass:: pyocto.associator.VelocityModel1D
    :members:
    :undoc-members:
    :show-inheritance:
