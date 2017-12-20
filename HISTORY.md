# Revision History for `matdb`

## Revision 0.0.7

- Renamed the databases to be groups and the `matdb/database/basic.py`
  `Database` class to be `Group`.
- Renamed the classes in `matdb/database/controller.py`,
  changed`DatabaseSequence` to just `Sequence` and `SequenceRepeater`
  to `Repeater`.
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
