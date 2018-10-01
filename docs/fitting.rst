Potential Fitting API
======================

At the moment, `matdb` only supports fitting using MTP. However, other
fitting routines may be added later.

YML File
--------

The fitting section of the yml file controls how the desired
interatomic potentials will be trained. The section contains the
fields `dbs`, `execution`, and `fits`.

- **dbs**: The `dbs` section of the yml file details the databases to
  be used for fitting. This can be over written later on in the `fits`
  section for a specific interatomic potential.
- **execution**: The `execution` section contains information needed
  for the job file and job submission.

fits
****

The `fits` section of the yml file contains details pertaining to each
individual fit desired. At present `matdb` only supports one type of
fit, howover, it has been built so that multiple fits can be made from
the same database so that different fitting methods can be
compared. The entries to this field are:

- **name**: The name of the fit.
- **dbs**: To overwrite the global `dbs` flag if desired.
- **steps**: Here the details for the fitting method are entered along
  with the actual fitting method to use via the `type` option.

Example yml File
****************

An example fitting section of the yml file.

.. code-block:: yaml
		
  fitting: 
     dbs: ['*'] 
     execution:
       template: 'run_mtp_ml.sh'
       time: 10
       ntasks: 1
       nodes: 1
       mem_per_cpu: 500MB
       job_name: 'AgPt-fitting'
     fits:
       - name: "AgPd"
	 dbs: ["*"]
	 steps:
           - type: "mtp.MTP"
             select: 
               selection-limit: 200
             species:
	       - "Ag"
	       - "Pd"

.. toctree::
   :maxdepth: 1

   fitting/mtp.rst
