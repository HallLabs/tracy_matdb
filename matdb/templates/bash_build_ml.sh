#!/bin/bash

# Iterator and run the jobs under all sub-directories
workdir=$(echo {{execution_path}} | awk -F'/[^/]*$' '{print $1}' )

for f in *; do
    if [ -d ${workdir}/${f} ]; then
        # ignore the configuration which has already been calculated
        # this is a temp fix for DFT_QE
        if [ -f ${workdir}/${f}/espresso.out ]; then
           continue
        fi

        cd ${workdir}/${f}
        # QE calculation will be timed out if runs up to the time out parameter in the yml file
        # Note that to turn off timeout, specify 0 for the exec_time_out_minutes parameter
        timeout -sTERM {{ exec_time_out_minutes }}m {{ exec_path }} -i espresso.pwi > espresso.out
    fi
done
