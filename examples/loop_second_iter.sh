#!/bin/bash

counter=2
while [ $counter -le 20 ]
do
  ./mtp_train_second_iter.sh
  cp ./compute/MTP/CoWV/CoWV/CoWV/mtp/train.cfg ./compute/MTP/CoWV/CoWV/CoWV/mtp/train.cfg_iter_$counter
  ((counter++))
  echo "pause for 5 minutes. we can safely kill the process during this time."
  sleep 5m
done

