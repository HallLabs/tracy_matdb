# Revision History for `matdb`


## Revision 1.1.6
- Working distortion groups.
- Added `matdb_find.py` script for testing patterns against controllers.
- Added QE calculator with 100% test coverage.
- Added the `stress_name` function to the basic calculator class.

## Revision 1.1.5
- Updated the Manual Group so that it can take an optional extractable
  flag which disapbles calulations when set to false.
- Updated `database/__init__.py' so that the user can input a
  transformation dictionary that will change the seed files passed in
  according to a transformation function.
- Created `can_extract()` method for the `Manual` group.
- Fixed missing `msg` import in `database/utility.py`.

## Revision 1.1.4
- Removed assignments to empty dicts from the function definitions as
  mention in Issue [#56](/../../issues/56).
- Added a warning when a database group being requested hasn't been
  implemented yet.
- Implemented the Prototypes database group.
- Added an exception for when the VASP executable can't be found when
  trying to determine the VASP version.

## Revision 1.1.3
- Changed the VASP calculator to use a series of POTCARs stored in the
  POTCARS dir.

## Revision 1.1.2
- Fixed the MTP method so that it correctly constructs the train.cfg
  file.
- Added the `Simple` database to the repo for when configurations need
  to be calculated without additional configurations being generated.
- Changed the VASP calculator to use a local instead of a global POTCAR.
- Fixed MTP method so that the POSCAR's get the zeros put back in
  after a calculation.
- Added new slurm template.
- Fixed more issues in the mtp commands.
- Implemented a re-write of the POTCAR by the mtp method that the ASE
  read can get the species from the title.
- Debugged Active Group setup.

## Revision 1.1.1

- Added the controll directory to the calling interface for the
  creation of calculators in the `Atoms` object.
- Added an check for `SinglePointCalculator` when reading the atoms
  from hdf5 files (we don't over-write the `ase`
  `SinglePointCalculator` so it dosen't have a `to_dict` method that
  would save the required fields so instead we just skip the
  calculator setup).
- Removed the folder argument from the calculators `to_dict` method
  since it isn't needed/used.
- Fixed `datetime` import in `matdb/database/__init__.py`.
- Removed duplicate `ran_seed` assignment in the database
  `Controller`.
- Fixed the Vasp calculator so that the `environ` variable get set in
  the `__init__` function.
- Fixed bugs in the atoms reading and writing routines.
- Added the Vasp version to the vasp stub.
- Moved `save_dict_to_h5` from matdb/utility.py to matdb/io.py.
- Fixed minor typos and bugs in the Vasp calculator.
- Removed duplicate `todict` method from the `Quip` calculotar.
- Added entry to kwargs dict of `to_dict` method of the `aflux`
  calculator.
- Fixed some minor typos in the `database` classes.
- Fixed the last bugs in the required packages script, include
  removing 'argparse' and 'setuptools' from the package list since
  'pip freeze' dosen't list them.
- Added an exception that gets printed if the 'Quip' calculator can't
  be loaded. This try and except fixes import errors if the caculator
  subpackage.
- Fixed bugs in the enumerated database introduced by updates to
  `phenum`.
- Fixed hashing problem in enumerated database.
- Added `pre_comp_atoms.h5` to creation of folders in the group
  class. It's removed after extraction has been completed.

## Revision 1.1.0

- The database finalize method has been implemented as described in
  Issue [#31](/../../issues/31).
- Rename the cleanup methods to extract and created a new cleanup
  method on the calculators that performs the desired level of cleanup
  as described in Issue [#35](/../../issues/35).

## Revision 1.0.9

- Fixed the bug reported in Issue [#33](/../../issues/33).
- Fixed atoms object after format change to calculators.
- Added progress bars to the database setup and cleanup methods as
  suggested in Issue [#26](/../../issues/26).
- Added the hash methods to group, database, and controller as well as
  the verify_hash method to the controller as described in Issue
  [#30](/../../issues/30).
- Added the time stamp to the second line of the uuid files. Also
  fixed some errors with the overwriting of objects in the group
  settings in which the new uuid for the new objects wasn't getting
  saved to file. This resolves Issue [#37](/../../issues/37).
- Implemented the `to_dict` method on the calculators as described in
  Issue [#29](/../../issues/29).
- Added the `contr_dir` to store tho controller directory for the
  calculators. This will be needed to implement Issue
  [#21](/../../issues/21).
- Implemented the creation of the `POTCAR` like file in the
  calculators instead of the `controller` as described in Issue
  [#21](/../../issues/21).
- Implemented the bug fix described in Issue [#22](/../../issues/22)
  so that `xc` only needs to be set in one place for the `vasp`
  calculator.
- Implemented dependency version checking and storing on matdb.

## Revision 1.0.8

- Moved all database specific functions from `matdb.utility` to
  `matdb.database.utility`.
- Removed the `split` method for the `legacy` database and `Database`
  classes and made a universal `split` function in
  `matdb.database.utility`.
- Moved the classes contained in `matdb.database.controller.py` and
  `matdb.database.basic.py` to `matdb.database.__init__.py`, and fixed
  relative imports.
- Switched random seeds to be universal when defined in controller.
- All calculators now take a random seed.

## Revision 1.0.7

- Added CFG support to `matdb.io` and increased unit test coverage.
- Adjusted dependencies list for `setup.py`.

## Major Revision 1.0

Includes a major, API-breaking refactor of the database handling for
`matdb`. We refactored the generation to include a recursive group
structure for handling multiple seeds, chained groups, and parameter
grids. Some features of the major revision are:

- Recursive group structure allows for parameter grids.
- Chaining of groups via a new `rset` property.
- Groups can define a set of training configurations.
- Scientific reproducibility enforced by having deterministic outcomes
  for a fixed set of input parameters to groups.
- `uuid`-based identifiers for all databases, groups and atoms
  objects.
- HDF5-based serialization for storage of configurations and
  calculator results.
- Hybrid support for ASE and quippy `Atoms` objects.
- Asynchronous calculator support for basic ASE calculators.

## Revision 0.0.7

- Implented the functions needed to indentify the seed configurations
  for seeded databases.
- Implemented the full Enumerated database for the new workflow.
- Added new functions to the `matdb/database/basic.py` `Group` class
  that allow the parameters passed to the function as part of a
  parameter grid to be saved to a json file in the folder in which the
  group with those parameters will be executed. Also added a function
  that allows the saved parameters to be read back in.
- Enabled looping over the parameters in the Parameter Grid in
  `matdb/database/basic.py`.
- Renamed the databases to be groups and the `matdb/database/basic.py`
  `Database` class to be `Group`.
- Renamed the classes in `matdb/database/controller.py`,
  changed`DatabaseSequence` to just `Sequence` and removed `SequenceRepeater`.
- Added the new subroutiens (`flatten_dict`,`special_values`,`special_functions`,
  `slicer`, `is_number` and `_py_exectue`) to `matdb/utilit.py` that are needed
  for the refactor.
- Implemented the ParametrGrid class in `matdb/database/controller.py`.
- Added struct_enum.out files for the gss search for `fcc`, `bcc`,
  `sc`, and `hcp` structures.

## Revision 0.0.6

- Replaced self.collections with self.seeded in `matdb/database/controller.py`.
- Assigned self.steps to the Sequence database.
- Added `support/matdb_sbatch.py` for the unit testing stub for sbatch.
- Added `matdb/querries.py` so that the different environments, i.e.,
  sbatch and Tracy, can be handled appropriately.
- Set the `self._bands` attribute in the phonon database.
- Added `EnumDatabase` for generating enumerated runs.
- Added `lattice.in` template for phenum's use.
- Fixed the `fileformat` in `matdb/kpoints.py` to be `vasp-ase` to
  match materialscloud server update.
- Changed `matdb/plotting/potentials.py` so that if the
  `atoms.config_type` attribute does not exist it doesn't throw
  errors.

## Revision 0.0.5

- Added `DynamicDatabase` for generating MD runs.

## Revision 0.0.4

- Debugged modulation of phonon database (using specified amplitudes).
- The inference of amplitudes is still not coded (is present in a prototyping notebook).

## Revision 0.0.3

- Debugged the modulation sub-config generation.
- Debugged the phonon amplitude calibration.
- Added status messages for database DFT execution progress.
- Debugged job array submission.
- Added several new statements to prevent re-building of setup and cleanup phases in
  databases that are already finished.

## Revision 0.0.2

- Added support for job file submission.
- Added execution and cleanup to the builder script interface (not debugged yet).

## Revision 0.0.1

- Added the script for interacting with database controller.
