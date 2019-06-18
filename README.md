# Tracy MatDB (Materials Database)

Code for generating configuration databases to train interatomic potentials.

The interatomic potentials are created using [MLIP](http://gitlab.skoltech.ru/shapeev/mlip). Which uses [Quantum Espresso](https://gitlab.com/QEF/q-e).

# Table of Contents

We describe step to run the project either on a local machine or a remote machine on AWS.

- [Introduction](#Introduction)
- [Build Containers Locally](#Build-Containers-Locally)
- [Running the Project](#Running-the-Project)
- [Remote Run Using AWS](#Remote-Run-Using-AWS)
    1. [Prepare AWS Instance](#Prepare-AWS-Instance)
    2. [Prepare and Upload Images](#Prepare-and-Upload-Images)
    3. [Setup Remote Machine](#Setup-Remote-Machine)
- [Miscellaneous](#Miscellaneous)
    - [API Documentation Generation](#API-Documentation-Generation)
        1. [Dependencies for Generation](#Dependencies-for-Generation)
        2. [Build API Documentation](#Build-API-Documentation)
        3. [Extract from Container](#Extract-from-Container)
    - [Unit Tests](#Unit-Tests)
    - [Helpful Docker Commands](#Helpful-Docker-Commands)
    - [Monitoring Tips](#Monitoring-Tips)
    - [Contributing](#Contributing)

# Introduction

# Build Containers Locally

<!-- Change URLs when project ownership changes. -->

Clone all the related projects:
- [Tracy Docker](https://github.com/NorimaConsulting/tracy_docker)
- [Tracy Matdb](https://github.com/NorimaConsulting/tracy_matdb)

Use `ssh` or `https` at your conviece. It is *recomended* that they be in a folder that makes it easy for you to find. It is **required** that the cloned projects be siblings of eachother.
- If you plan to use `ssh` keys, there is a good guide from [github](https://help.github.com/en/articles/adding-a-new-ssh-key-to-your-github-account).

The resulting folder structure should look like:

```
tracy_root
├── tracy_docker
└── tracy_matdb
```

Follow the instructions in the [Tracy Docker readme](https://github.com/NorimaConsulting/tracy_docker/blob/master/README.md) to build all relevent prerequisite docker images.

<!-- matdb               latest   ...          ...          1.87GB -->

To verify, run `docker images`. You should see similar output:
```
REPOSITORY          TAG      IMAGE ID     CREATED      SIZE
mlip4tracy          latest   ...          ...          1.98GB
pslibrary4tracy     latest   ...          ...          6.56GB
dft_qe4tracy        latest   ...          ...          2.43GB
ubuntu4tracy        16.04    ...          ...          1.54GB
ubuntu              16.04    ...          ...          119MB
```

You may see more images, but the above should be present.

We can now build the `matdb` image. Using a shell of your choice go into the `tracy_matdb` directory and run the `build_docker_image.sh` script to build the `matdb` image. Note that this image builds on top of `mlip4tracy`.

Now you can either run them locally, or on AWS.

# Running the Project

Some optional tips, you may want to use a terminal multiplexer like `tmux`. We have found it useful, especially if the host machine is remote.

Run the docker image with:
```bash
docker run -it --rm matdb /bin/bash
```

At this point you need to provide a `yml` file or use the example one. The format of the `yml` file is specified in the API Documentation (which must be first generated before viewing).

Then within the instance:
1. move the `/root/codes/` directory
2. run `./mtp_train.sh <yml-file-location>`
    - If you do not specify a `yml` file location then matdb will run based off the example `yml` file.
    - Note that the scripts run cell sizes of 3 or 4 runs until the `mtp` file has converged suitabliy then increases the cell size.

# Remote Run Using AWS

## Prepare AWS Instance

<!-- TODO: Incorperate deploy readme here -->

### generate an new ec2 instance
- instructions [here](./aws_deploy.md)
- a linux based one
- needs 64 gb of HD
- 8 gb of RAM
- create a key or use an existing

## Prepare and Upload Images

Note for this step, your local machine requires sufficient hard drive space.

```bash
➜ cd tracy_matdb/scripts/deploy
➜ ./deploy.sh # this extracts the the required images as jar files
```

Then you can upload the `.jar` files to the remote machine.

Within the `tracy_matdb/scripts/deploy` folder run the following command.
```bash
➜ IDENTITY_FILE="<somepath>/tracy-science-test-001.pem" \
    DESTINATION_URL="ec2-user@<ec2-id>.<region>.compute.amazonaws.com" \
    DESTINATION_PATH="~/" \
    ./deploy_aws.sh
```
Fill in the above information to facilitate a connection to your remote aws instance. If your `DESTINATION_URL` format differs, just use the one you are give from the `aws` console.

## Setup Remote Machine

First you need to make sure that `docker` is installed.

If it isn't, the command depends on the package manager on your system. The following are examples with `yum` and `apt`.
```bash
➜ sudo yum -y install docker
```
```bash
➜ sudo apt-get install docker
```

Next you need to make sure that the `docker-service` is running.

```bash
➜ sudo service docker start
```

Then load the docker images into docker service. The files and script should be in the `$HOME` directory if following previous steps.
```bash
➜ cd $HOME
➜ ./load_docker_images_aws.sh
```

After the images are loaded, you can run the project following the steps in the [run section](#Running-the-Project).

# Miscellaneous

## Build Documentation

### Dependencies for Generation
- `python3 -m pip install git+https://github.com/sphinx-doc/sphinx`
- `python3 -m pip install sphinxcontrib-napoleon`

<!-- Clarify steps, which are in the container, which are on host machine -->

### API Documentation Generation

Move to the `tracy_matdb/docs` direcotry. Then run `make` to see what format you want to make the documentation as.
    - To create an `html` version run `make html`.

The following is an example:
```bash
➜ cd tracy_matdb/docs
➜ make html
```

### Extract from Container

You want to extract the relevent files in the `_build` directory that was recently generated.

You can use `docker cp` to extract them from inside a docker image.
```bash
➜ docker cp <container-id>:/root/codes/matdb/docs/_build <host-destination>
```
- the `<container-id>` can be found using `docker ps`.
- `<host-destination>` is a file location on your local machine that you want to place the files.

The extracted files can then be viewed. If you created an `html` version, you can open `<host-destination>/_build/html/index.html` in your web broswer of choice.

## Unit Tests
From within an instance of the docker image `matdb` run the following:
```bash
➜ python3 -m pytest /root/codes/matdb/tests
```

To run the tests with a code coverage tool.
- Install `pytest-cov` with `python3 -m pip install pytest-cov`.
Then to run the unit tests with the coverage tool:
```bash
➜ python3 -m pytest --cov=/root/codes/matdb /root/codes/matdb/tests
```

## Helpful Docker Commands

Depending on how `docker` is installed, you made need to run commands with `sudo`.

Stop all running container
```bash
docker stop $(sudo docker ps -a -q)
```

Remove all docker images
```bash
docker rmi $(docker image ls -a  --format  "{{.ID}}")
```

Remove stop all containers and remove unused images
```bash
docker system prune
```

Create a second connection to an existing container
```bash
➜ docker exec -it <container-id> /bin/bash
```
The container ID can be found using `docker ps`

## Monitoring Tips

## Intermediate Files

- `to-relax.cfg`
    - This is the file contains structures needed to be relaxed. The idea is that the IAP should be able to relax all the structures listed in this file. If a structure can not be relaxed, it then is added to the `new_training.cfg` file and is ultimately added to the training set(`train.cfg`).
    - This file is generated at the first iteration for each atom cell iteration.
- `new_training.cfg(new_training.cfg_iter_?)`
    - Each iteration will generate some new structures which couldn’t be relaxed by the current IAP.  These new structures will be added to the training set (train.cfg)  at the beginning of the next iteration. I saved a copy of this file for each iteration for debug purpose. For example: new_training.cfg_iter_6 is for the 6th iteration. If this file is empty, it means it converges at the iteration this file is corresponding to. For example, if `new_training.cfg_iter_17` is empty, it shows it converged at iteration 17. Notice that, as long as the number of new structures generated at an iteration is less than or equal to the next_cell_threshold defined in the `yml` file, it considered converged
- `train.cfg(train.cfg_iter_?)`
    - This is the training set.  I saved a copy of this file for each iteration for debug purpose.
- `pot.mtp(pot.mtp_iter_?)`
    - This is the potentials.  I saved a copy for each iteration for debug purpose.
- `training.txt (training.txt_iter_?)`
    - This is the log file for the `mtp train` process.  At the bottom of the file, it shows the training errors which Wiley would be interested in. Especially for the `Energy per atom`.
- `status.txt`: this contains status code for each step in an iteration.
    - Some of the status:
        ```python
        "relax_setup {0} {1}".format(self.iter_count, self.cell_iter)
        "relax {0} {1}".format(self.iter_count, self.cell_iter)
        "select {0} {1}".format(self.iter_count, self.cell_iter)
        "add {0} {1}".format(self.iter_count, self.cell_iter)
        "done {0} {1} {2}".format(self.iter_count, self.cell_iter, len(new_configs))
        ```
        - Refer to `command()` method in `fitting.mtp.py module` for a complete status and it’s meaning.
- `jobfile.sh`
    - This file contains `mtp` command to be executed. Each iteration has it’s own `mtp` commands need to be carried out.
- `iter_?.pkl` files in `Active` database.
    - The `Active` database resides at `/root/codes/compute/MTP/CoWV/Active/active.CoWV` for our example `CoWV` structures.  Each iteration will have it’s `pkl` file generated at the `Active` database root directory. Each `pkl` file contains the new structures for the specific iteration. Notice that, the number of structures in `new_training.cfg_iter_?` and `iter_?.pkl` should be the same. For example, `new_training.cfg_iter_9` and `iter_9.pkl` should have the same number of structures. But `train.cfg_iter_10` minus `train.cfg_iter_9` might have less structures then in the two files. That is because QE calculation could fail on some of the new structures.
- `matdb/templates/bash_build_ml.sh`
    - the template file used to generate the `jobfile.sh` in the `Active` database root directory defined above.  This template file works only for `QE` calculation.
- `pkl` means is the extension for a python pickle file.


## Contributing

Before contributing to `matdb` please read the [contributing guidelines](/CONTRIBUTING.md).
