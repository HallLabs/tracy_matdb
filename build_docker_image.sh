#! /bin/bash

DEV_MODE=false
while [[ $# -gt 0 ]] && [[ ."$1" = .--* ]] ;
do
    opt="$1";
    shift;              #expose next argument
    case "$opt" in
        "--devmode" )
           DEV_MODE=true;;
        *) echo >&2 "Invalid option: $@"; exit 1;;
   esac
done

# build matdb
# the purpose of doing export/import is to get rid of the annoying waring "Unexpected end of /proc/mounts line"
if [ ${DEV_MODE} == true ] ; then
  docker image build -t matdb_stable . --build-arg DEV_MODE="YES"
  if [ $? -ne 0 ]; then
    echo "ERROR on building dev_mode matdb"
    exit 1
  fi
else
  docker image build -t matdb_stable . --build-arg DEV_MODE="NO"
  if [ $? -ne 0 ]; then
    echo "ERROR on building matdb"
    exit 1
  fi
fi

docker run -d --name=matdb_temp matdb_stable
if [ $? -ne 0 ]; then
  echo "ERROR on running matdb"
  exit 1
fi

docker export matdb_temp | docker import - matdb_stable
if [ $? -ne 0 ]; then
  echo "ERROR on exporting/importing matdb_temp"
  exit 1
fi

docker stop matdb_temp
if [ $? -ne 0 ]; then
  echo "ERROR on stopping matdb_temp"
  exit 1
fi

docker rm matdb_temp
if [ $? -ne 0 ]; then
  echo "ERROR on removing container matdb_temp"
  exit 1
fi

# run the docker image as a service, enabled gdb from docker 
docker run --cap-add=SYS_PTRACE --security-opt seccomp=unconfined -it --rm -d matdb_stable /bin/bash
