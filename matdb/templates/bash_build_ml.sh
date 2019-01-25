#!/bin/bash

# Iterator and run the jobs under all sub-directories
workdir=$(echo {{execution_path}} | awk -F'/[^/]*$' '{print $1}' )

for f in *; do
    if [ -d ${workdir}/${f} ]; then
        cd ${workdir}/${f}
        {{ exec_path }} -i espresso.pwi > espresso.out
    fi
done
