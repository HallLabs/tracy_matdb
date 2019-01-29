#!/bin/bash

# Iterator and run the jobs under all sub-directories
workdir=$(echo {{execution_path}} | awk -F'/[^/]*$' '{print $1}' )

for f in *; do
    if [ -d ${workdir}/${f} ]; then
        # ignore the configuration which has already been calculated
        # this is temp fix for DFT_QE
        if [ -f ${workdir}/${f}/espresso.out ]; then
           continue
        fi

        cd ${workdir}/${f}
        {{ exec_path }} -i espresso.pwi > espresso.out
    fi
done
