## Quick Start Instructions

### 1 - generate an new ec2 instance
- a linux based one
- needs 64 gb of HD
- 8 gb of RAM
- create a key or use an existing
- instructions [here](./aws_deploy.md)

### 2 - recomended setup of a github
- https://help.github.com/en/articles/adding-a-new-ssh-key-to-your-github-account

clone all the projects
https://github.com/NorimaConsulting/tracy_matdb
https://github.com/NorimaConsulting/tracy_docker
https://github.com/NorimaConsulting/tracy_science

**NOTE** change to docker-setup branch on all of them. (will be removed later)

### 3 - build everything (in this order)
```bash
➜ cd tracy_docker/ubuntu16.04/ && ./build_docker_images.sh
➜ cd tracy_science/mtp && ./build_docker_image.sh
➜ cd tracy_matdb && ./build_docker_image.sh
```

### 4 - get everything setup for deployment
```bash
➜ cd /tracy_matdb/scripts/deploy
➜ run deploy.sh
➜ IDENTITY_FILE="<somepath>/tracy-science-test-001.pem" DESTINATION_URL="ec2-user@<ec2-id>.<region>.compute.amazonaws.com" DESTINATION_PATH="~/" ./deploy_aws.sh
```

### 5 - connect to remote machine
```bash
➜ ssh -i <identity-file> ec2-user@<ec2-id>.<region>.compute.amazonaws.com
```

### 6 - install docker on remote machine
```bash
➜ sudo yum -y install docker
➜ sudo apt-get install docker
```

### 7 - run the docker service
```bash
➜ sudo service docker start
```

### 8 - load the docker images into the docker service
```bash
➜ load_docker_images_aws.sh
```

### 9 Optional - install a terminal multiplexer (if you want to run things over multiple days)
```bash
➜ yum install tmux
```

#### 9.1 - with tmux, create a new session
```bash
➜ tmux new
```

#### 9.2 - attach to it later with
```bash
➜ tmux a <session-index>
```

### 10 - run the docker image
```bash
➜ sudo docker run -it --rm matdb /bin/bash
➜ cd /root/codes/
➜ ./mtp_train.sh
```

### 11 - Optional: With a different connection, you can connect to the same container
```bash
➜ sudo docker ps
➜ sudo docker exec -it e55cb9787958 /bin/bash
```
