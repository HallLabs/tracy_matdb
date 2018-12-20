#!/bin/bash

~/mtp_train_first_iter.sh
if [ $? -ne 0 ]; then
  echo "ERROR on running first iteration"
  exit 1
fi

echo "pause for 5 minutes. we can safely kill the process during this time."
sleep 5m

counter=2
while [ $counter -le 20 ]
do
  ~/mtp_train_second_iter.sh
  if [ $? -ne 0 ]; then
    echo "ERROR on running iteration $counter"
    exit 1
  fi

  cp ~/compute/MTP/CoWV/CoWV/CoWV/mtp/train.cfg ~/compute/MTP/CoWV/CoWV/CoWV/mtp/train.cfg_iter_$counter
  ((counter++))
  echo "pause for 5 minutes. we can safely kill the process during this time."
  sleep 5m
done

