Sum of creep laws
=================

Overview
--------

This object provides the sum of several individual scalar creep laws

.. math::
   \dot{\varepsilon}^{cr} = \sum_{i=1}^{n} \dot{\varepsilon}^{cr}_{i} 

where each of :math:`\dot{\varepsilon}^{cr}_{i}` represent individual scalar
creep laws.

Parameters
----------

.. csv-table::
   :header: "Parameter", "Object type", "Description", "Default"
   :widths: 12, 30, 50, 8

   ``models``, :cpp:class:`std::vector<neml::ScalarCreepRule>`, Vector of creep models, No

Class Description
-----------------

.. doxygenclass:: neml::SumScalarCreep
   :members:
   :undoc-members:
