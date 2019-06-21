#!/bin/bash

##
## This script assumes that we have only one fitter configured in the yaml specification file in the fitting.fits section
##
## You may specify the YAML specification file from the command line as a parameter like this:
##     mtp_train.sh Qe_CoWV.yml
## Or it will find one from the current folder if you don't specify one.
##

## Function to extract values from the YAML specification file
yaml() {
    python3 -c "exec(\"import yaml\ntry:\n\tprint(yaml.load(open('$1'))$2)\nexcept:\n\tprint("0")\")"
}

## Function to save some intermediate output files for log/debug purpose
saveOutputFiles() {
    cp ${ROOT_DIR}/${FIT_NAME}/${FIT_NAME}/mtp/train.cfg ${ROOT_DIR}/${FIT_NAME}/${FIT_NAME}/mtp/train.cfg_iter_$ITER
    cp ${ROOT_DIR}/${FIT_NAME}/${FIT_NAME}/mtp/pot.mtp ${ROOT_DIR}/${FIT_NAME}/${FIT_NAME}/mtp/pot.mtp_iter_$ITER
    cp ${ROOT_DIR}/${FIT_NAME}/${FIT_NAME}/mtp/training.txt ${ROOT_DIR}/${FIT_NAME}/${FIT_NAME}/mtp/training.txt_iter_$ITER
}

# Function to run the first iteration
firstIter() {
    ((ITER++))

    echo "Iteration $ITER starting at: `date --iso-8601=seconds`"
    ./mtp_train_first_iter.sh ${SPEC_FILE_BASE}
    if [ $? -ne 0 ]; then
      echo "ERROR on running first ITER"
      exit 1
    fi
    echo "Iteration $ITER ended at: `date --iso-8601=seconds`"

    # copy/save output files
    saveOutputFiles

    echo "pause for 5 minutes. we can safely kill the process during this time."
    sleep 5m
}

# Function to run the second iteration. The second itertion will be repeated until it converges
secondIter() {
    ((ITER++))

    echo "Iteration $ITER starting at: `date --iso-8601=seconds`"
    ./mtp_train_second_iter.sh ${SPEC_FILE_BASE}
    if [ $? -ne 0 ]; then
      echo "ERROR on running first ITER"
      exit 1
    fi
    echo "Iteration $ITER ended at: `date --iso-8601=seconds`"

    # copy/save output files
    saveOutputFiles

    echo "pause for 5 minutes. we can safely kill the process during this time."
    sleep 5m
}

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

# grab root directory from the yaml file
ROOT_DIR=$(yaml "$SPEC_FILE" "['root']")
eval ROOT_DIR=$ROOT_DIR

# grab the next_cell_threshold from the yaml file
NEXT_CELL_THRESHOLD=$(yaml "$SPEC_FILE" "['fitting']['fits'][0]['steps'][0]['next_cell_threshold']")

# grab the name of the fit
FIT_NAME=$(yaml "$SPEC_FILE" "['fitting']['fits'][0]['name']")

# build the initial database
./mtp_build.sh ${SPEC_FILE_BASE}

# initialize the current iteration
ITER=0

# for each cell iteration, run the first iteration once
firstIter

# repeat the second iteration until it's converged
while [ $ITER -le 50 ]
do
  database_size=`find "${ROOT_DIR}/${FIT_NAME}" -name "atoms.h5" | wc -l`
  training_size_before=`grep BEGIN ${ROOT_DIR}/${FIT_NAME}/${FIT_NAME}/mtp/train.cfg | wc -l`

  # if the database size and the training set size are not equal, report an error and exit
  if [ ${database_size} -ne ${training_size_before} ]; then
    echo "ERROR on populating training set."
    exit 1
  fi

  # run the second iteration
  secondIter

  new_configuration=`grep BEGIN ${ROOT_DIR}/${FIT_NAME}/${FIT_NAME}/mtp/new_training.cfg_iter_$ITER | wc -l`

  # if the new configuration is less than or equal to the next cell threshold, move on to the next cell iteration.
  if [ ${new_configuration} -le ${NEXT_CELL_THRESHOLD} ]; then
    firstIter
  fi

done
