#!/bin/bash

COUNTER=1
STEP=1
echo "Step $STEP starting at: `date --iso-8601=seconds`"
matdb_train.py Qe_CoWV -t
if [ $? -ne 0 ]; then
  echo "mtp_train_first_iteration step $STEP Errored"
  exit 1
fi

while [ $COUNTER -lt 5 ]; do
  let STEP+=1
  echo "Step $STEP starting at: `date --iso-8601=seconds`"
  matdb_train.py Qe_CoWV -x
  if [ $? -ne 0 ]; then
    echo "mtp_train_first_iteration step $STEP Errored"
    exit 1
  fi

  let STEP+=1
  echo "Step $STEP starting at: `date --iso-8601=seconds`"
  matdb_train.py Qe_CoWV -t
  if [ $? -ne 0 ]; then
    echo "mtp_train_first_iteration step $STEP Errored"
    exit 1
  fi

  let COUNTER+=1
done

