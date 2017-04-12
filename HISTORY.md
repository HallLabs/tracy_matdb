# Revision History for `matdb`

## Revision 0.0.4

- Debugged modulation of phonon database (using specified amplitudes).
- The inferance of amplitudes is still not coded (is present in a prototyping notebook).

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