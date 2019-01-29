#!/bin/bash

COUNTER=1
STEP=1
echo "Step $STEP starting at: `date --iso-8601=seconds`"
matdb_build.py Qe_CoWV -e --verbose 5
if [ $? -ne 0 ]; then
  echo "mtp_train_second_iteration step $STEP Errored"
  exit 1
fi

while [ $COUNTER -lt 4 ]; do
  let STEP+=1
  echo "Step $STEP starting at: `date --iso-8601=seconds`"
  matdb_train.py Qe_CoWV -t --verbose 5
  if [ $? -ne 0 ]; then
    echo "mtp_train_second_iteration step $STEP Errored"
    exit 1
  fi

  let STEP+=1
  echo "Step $STEP starting at: `date --iso-8601=seconds`"
  matdb_train.py Qe_CoWV -x --verbose 5
  if [ $? -ne 0 ]; then
    echo "mtp_train_second_iteration step $STEP Errored"
    exit 1
  fi

  let COUNTER+=1
done

let STEP+=1
echo "Step $STEP starting at: `date --iso-8601=seconds`"
matdb_train.py Qe_CoWV -t --verbose 5
if [ $? -ne 0 ]; then
  echo "mtp_train_second_iteration step $STEP Errored"
  exit 1
fi
