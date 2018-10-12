#! /bin/bash

# build matdb
# the purpose of doing export/import is to get rid of the annoying waring "Unexpected end of /proc/mounts line"
docker image build -t matdb .
docker run -d --name=matdb_temp matdb
docker export matdb_temp | docker import - matdb
docker stop matdb_temp
docker rm matdb_temp

# run the docker image, enable gdb from docker 
docker run --cap-add=SYS_PTRACE --security-opt seccomp=unconfined -it --rm matdb /bin/bash

