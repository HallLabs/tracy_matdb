#!/bin/bash

echo "first iteration step 1"
matdb_train.py Qe_CoWV -t 
if [ $? -ne 0 ]; then
  echo "mtp_train_first_iteration step 1 Errored"
  exit 1
fi

COUNTER=1
STEP=1
while [ $COUNTER -lt 5 ]; do
  let STEP+=1
  echo "first iteration step $STEP"
  matdb_train.py Qe_CoWV -x 
  if [ $? -ne 0 ]; then
    echo "mtp_train_first_iteration step $STEP Errored"
    exit 1
  fi

  let STEP+=1
  echo "first iteration step $STEP" 
  matdb_train.py Qe_CoWV -t 
  if [ $? -ne 0 ]; then
    echo "mtp_train_first_iteration step $STEP Errored"
    exit 1
  fi

  let COUNTER+=1
done

