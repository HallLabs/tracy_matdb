#!/bin/bash

workdir="{{execution_path}}"
workdir=${workdir::-2}

# Iterator and run the jobs under E.* directories 
for dir in E.*
do
  cd "${workdir}${dir}" 
  {{ exec_path }} -i {{ input_file }} > {{ output_file }}
done
