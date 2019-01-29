#!/bin/bash

matdb_build.py Qe_CoWV -s --verbose 5
if [ $? -ne 0 ]; then
  echo "mtp_build step 1 Error"
  exit 1
fi

matdb_build.py Qe_CoWV -x --verbose 5
if [ $? -ne 0 ]; then
  echo "mtp_build step 2 Error"
  exit 1
fi

matdb_build.py Qe_CoWV -e --verbose 5
if [ $? -ne 0 ]; then
  echo "mtp_build step 3 Error"
  exit 1
fi
