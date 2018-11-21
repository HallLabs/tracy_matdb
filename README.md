

# Materials Database

Code for generating configuration databases to train interatomic potentials. [API Documentation](https://rosenbrockc.github.io/matdb/)

## Contributing

Before contributing to `matdb` please read the Contributing
[guidelines](https://github.com/rosenbrockc/matdb/blob/master/CONTRIBUTING.md).

## Docker for Tracy

* Download code for `docker` project, run `build_docker_images.sh` script under `docker/ubuntu16.04` folder to build base iamges
* Download code for `mtp` project, run `build_docker_image.sh` script under `mtp` folder to build the mtp image
* Download code for `matdb` project, run `build_docker_image.sh` script under `matdb` folder to build the matdb image
* Can specify option "--devmode" to keep source code in the docker image when running the build scripts.
