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

COUNTER=1
STEP=1
echo "Step $STEP starting at: `date --iso-8601=seconds`"
matdb_train.py ${SPEC_FILE_BASE} -t --verbose 5
if [ $? -ne 0 ]; then
  echo "mtp_train_first_iteration step $STEP Errored"
  exit 1
fi

while [ $COUNTER -lt 5 ]; do
  let STEP+=1
  echo "Step $STEP starting at: `date --iso-8601=seconds`"
  matdb_train.py ${SPEC_FILE_BASE} -x --verbose 5
  if [ $? -ne 0 ]; then
    echo "mtp_train_first_iteration step $STEP Errored"
    exit 1
  fi

  let STEP+=1
  echo "Step $STEP starting at: `date --iso-8601=seconds`"
  matdb_train.py ${SPEC_FILE_BASE} -t --verbose 5
  if [ $? -ne 0 ]; then
    echo "mtp_train_first_iteration step $STEP Errored"
    exit 1
  fi

  let COUNTER+=1
done

