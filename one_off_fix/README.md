

# This folder contains temporary fixes for third-party libraries. 

## Fixes

* pip 
This is to fix the annoying false error messages when the matdb trying to check for depends.

* espresso.py 
This is to fix an issue in ASE library.
ASE library assumes each attribute of PP_HEADER element in the ps file is in a separated line
but pslibrary generates all the attributes in a single line. 
See line 1138-1139, in method grep_valence()
