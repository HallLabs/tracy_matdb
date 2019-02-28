#!/bin/bash

# grab root directory from the yml file
ROOT_DIR=`awk '/'root'/ {print $2}' *.yml | tr -d "\'" `
eval ROOT_DIR=$ROOT_DIR

# build the initial database
./mtp_build.sh

# for each cell iteration, run the first iteration once
ITER=1
echo "Iteration $ITER starting at: `date --iso-8601=seconds`"
./mtp_train_first_iter.sh
if [ $? -ne 0 ]; then
  echo "ERROR on running first ITER"
  exit 1
fi
echo "Iteration $ITER ended at: `date --iso-8601=seconds`"

# copy/save output files
cp ${ROOT_DIR}/CoWV/CoWV/mtp/train.cfg ${ROOT_DIR}/CoWV/CoWV/mtp/train.cfg_iter_$ITER
cp ${ROOT_DIR}/CoWV/CoWV/mtp/pot.mtp ${ROOT_DIR}/CoWV/CoWV/mtp/pot.mtp_iter_$ITER
cp ${ROOT_DIR}/CoWV/CoWV/mtp/training.txt ${ROOT_DIR}/CoWV/CoWV/mtp/training.txt_iter_$ITER

echo "pause for 5 minutes. we can safely kill the process during this time."
sleep 5m

# repeat the second iteration until it's converged
while [ $ITER -le 50 ]
do
  database_size=`find compute -name "atoms.h5" | wc -l`
  training_size_before=`grep BEGIN ${ROOT_DIR}/CoWV/CoWV/mtp/train.cfg | wc -l`

  # if the database size and the training set size are not equal, report an error and exit
  if [ ${database_size} -ne ${training_size_before} ]; then
    echo "ERROR on populating training set."
    exit 1
  fi

  ((ITER++))
  echo "Iteration $ITER starting at: `date --iso-8601=seconds`"
  ./mtp_train_second_iter.sh
  if [ $? -ne 0 ]; then
    echo "ERROR on running ITER $iteration"
    exit 1
  fi
  echo "Iteration $ITER ended at: `date --iso-8601=seconds`"

  # copy/save output files
  cp ${ROOT_DIR}/CoWV/CoWV/mtp/train.cfg ${ROOT_DIR}/CoWV/CoWV/mtp/train.cfg_iter_$ITER
  cp ${ROOT_DIR}/CoWV/CoWV/mtp/pot.mtp ${ROOT_DIR}/CoWV/CoWV/mtp/pot.mtp_iter_$ITER
  cp ${ROOT_DIR}/CoWV/CoWV/mtp/training.txt ${ROOT_DIR}/CoWV/CoWV/mtp/training.txt_iter_$ITER

  training_size_after=`grep BEGIN ${ROOT_DIR}/CoWV/CoWV/mtp/train.cfg | wc -l`

  # if the training size didn't change, assume it's converged. 
  if [ ${training_size_before} -eq ${training_size_after} ]; then
    echo "Converged."
    exit 0
  fi

  echo "pause for 5 minutes. we can safely kill the process during this time."
  sleep 5m
done

