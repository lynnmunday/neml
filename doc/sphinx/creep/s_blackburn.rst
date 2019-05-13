Blackburn creep law
===================

Overview
--------

This object implements a single Blackburn creep model with the form

.. math::
   \dot{\varepsilon}^{cr} = r \left(\varepsilon_{r} - \varepsilon^{cr} \right)

for temperature dependent parameters :math:`r` and :math:`\varepsilon_{r}`.

Parameters
----------

.. csv-table::
   :header: "Parameter", "Object type", "Description", "Default"
   :widths: 12, 30, 50, 8

   ``r``, :cpp:class:`neml::Interpolate`, Rate constant, No
   ``eps_r``, :cpp:class:`neml::Interpolate`, Saturation strain, No

Class description
-----------------

.. doxygenclass:: neml::BlackburnCreep
   :members:
   :undoc-members:
