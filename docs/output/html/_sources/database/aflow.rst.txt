AFLOW Database Documentation
============================

Example yml File Entry
**********************

.. code-block:: yaml

   steps:
     - type: "aflow.Aflow"
       catalog: 'icsd'
       batch_size: 100
       filters:
         - - - ["Egap", ">", 0]
	     - "&"
	     - ["Egap", "<", 2]
  	   - "|"
	   - - ["Egap", ">", 5]
	     - "&"
	     - ["Egap", "<", 7]
	 - - - ["species", "==" ,"Ag"]
	     - "&"
	     - ["species", "==", "Pd"]
       select:
         - "agl_thermal_conductivity_300K"
       orderby:
         keyword: "agl_thermal_conductivity_300K"
	 desc: True

AFLOW Database
**************

.. automodule:: matdb.database.aflux
   :synopsis: Database sub-classes for querying the `aflowlib.org`
              database.
   :members:      

