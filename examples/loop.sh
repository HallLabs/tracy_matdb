#!/bin/bash

ITER=1
echo "Iteration $ITER starting at: `date --iso-8601=seconds`"
~/mtp_train_first_iter.sh
if [ $? -ne 0 ]; then
  echo "ERROR on running first ITER"
  exit 1
fi
echo "Iteration $ITER ended at: `date --iso-8601=seconds`"

echo "pause for 5 minutes. we can safely kill the process during this time."
sleep 5m

while [ $ITER -le 20 ]
do
  ((ITER++))
  echo "Iteration $ITER starting at: `date --iso-8601=seconds`"
  ~/mtp_train_second_iter.sh
  if [ $? -ne 0 ]; then
    echo "ERROR on running ITER $iteration"
    exit 1
  fi
  echo "Iteration $ITER ended at: `date --iso-8601=seconds`"

  cp ~/compute/MTP/CoWV/CoWV/CoWV/mtp/train.cfg ~/compute/MTP/CoWV/CoWV/CoWV/mtp/train.cfg_iter_$ITER
  echo "pause for 5 minutes. we can safely kill the process during this time."
  sleep 5m
done

