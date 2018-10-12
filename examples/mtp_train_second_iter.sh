#!/bin/bash

echo "second iteration step 1"
matdb_build.py Qe_CoWV -e 
if [ $? -ne 0 ]; then
  echo "mtp_train_second_iteration step 1 Errored"
  exit 1
fi

COUNTER=1
STEP=1
while [ $COUNTER -lt 4 ]; do
  let STEP+=1
  echo "second iteration step $STEP"
  matdb_train.py Qe_CoWV -t
  if [ $? -ne 0 ]; then
    echo "mtp_train_second_iteration step $STEP Errored"
    exit 1
  fi

  let STEP+=1
  echo "second iteration step $STEP" 
  matdb_train.py Qe_CoWV -x 
  if [ $? -ne 0 ]; then
    echo "mtp_train_second_iteration step $STEP Errored"
    exit 1
  fi

  let COUNTER+=1
done

echo "second iteration step 8"
matdb_train.py Qe_CoWV -t
if [ $? -ne 0 ]; then
  echo "mtp_train_second_iteration step 8 Errored"
  exit 1
fi
