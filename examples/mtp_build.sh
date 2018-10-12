#!/bin/bash

matdb_build.py Qe_CoWV -s
if [ $? -ne 0 ]; then
  echo "mtp_build step 1 Error"
  exit 1
fi

matdb_build.py Qe_CoWV -x
if [ $? -ne 0 ]; then
  echo "mtp_build step 2 Error"
  exit 1
fi

matdb_build.py Qe_CoWV -e
if [ $? -ne 0 ]; then
  echo "mtp_build step 3 Error"
  exit 1
fi
