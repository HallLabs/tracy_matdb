

# Materials Database

Code for generating configuration databases to train interatomic potentials. [API Documentation](https://rosenbrockc.github.io/matdb/)

## Contributing

Before contributing to `matdb` please read the Contributing
[guidelines](https://github.com/rosenbrockc/matdb/blob/master/CONTRIBUTING.md).

## Build docker images for matdb 

* Download code from the `tracy_docker` project. This project contains docker files to build base images by `matdb`. Run `build_docker_images.sh` script under `ubuntu16.04` folder to build these base iamges
* Download code from the `mtp` project. It contains a docker file and a script file used to build the `mlip` docker image. Copy the two files to your `mlip` source code root directory. Then run `build_mlip_docker_image.sh` to build the `mlip` image
* Download code from the `tracy_matdb` project. Run `build_matdb_docker_image.sh` script to build the `matdb` docker image
* You may specify option "--devmode" to the build scripts to keep source code in the docker images
* Notice that we use `python3` to run `matdb`. `python 2.7` is also installed as `python` because it is required by third-party packages 

## Run matdb from a docker container
* After built the `matdb` docker image, you may start a docker container by typing `docker run -it --rm matdb /bin/bash`
* Then change to the directory `/root/codes` from inside the docker container
* You may change the `Qe_CoWV.yml` specification file according to your needs (Please refer to the API documentation mentioned above for making any changes to the `yml` file)
* Start the training process by running the `mtp_train.sh` script
* Wait patiently for the generating of the potential

## Unit tests
* From a docker container, change directory to `/root/codes/matdb`, type `python3 -m pytest .` to run a full unit test 
