#!/bin/bash

# grab the YAML specification file from the command line argument if it's specified
if [ $# -gt 0 ]; then
  SPEC_FILE=$1
fi

# grab the YAML specification file from the current directory if it's not specified from the command line
if [ ! -f "$SPEC_FILE" ]; then
  SPEC_FILE=$(ls *.yml 2>/dev/null)
fi

# exit if a YAML specification file could not be found
if [ ! -f "$SPEC_FILE" ]; then
  echo "No YAML specification file."
  exit 1
fi

# grab the base name of the YAML specification file
SPEC_FILE_BASE=$(echo "$SPEC_FILE" | cut -f 1 -d '.')

matdb_build.py ${SPEC_FILE_BASE} -s --verbose 5
if [ $? -ne 0 ]; then
  echo "mtp_build step 1 Error"
  exit 1
fi

matdb_build.py ${SPEC_FILE_BASE} -x --verbose 5
if [ $? -ne 0 ]; then
  echo "mtp_build step 2 Error"
  exit 1
fi

matdb_build.py ${SPEC_FILE_BASE} -e --verbose 5
if [ $? -ne 0 ]; then
  echo "mtp_build step 3 Error"
  exit 1
fi
