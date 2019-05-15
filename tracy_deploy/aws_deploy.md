# AWS Deployment Documentation:
# EC2 Instance Creation and How to Connect

1. Login to your AWS management console.
2. Tap EC2 to get into Resources screen, then select Running Instances, as show in this screenshot: [AWS Console](./aws-screenshot-1.png)
3. Now Tap Launch Instance, then follow the steps to create new AWS EC2 instance:[AWS Console - Create EC2 Instance](./aws-screenshot-3.png)
4. At the end of step "Review Instance Launch", create a new key pair, choose a key pair name, for example: matdb-aws-ec2, then tap "Download key Pair", you will have a pem file matdb-aws-ec2.pem stored to your local machine.[AWS Console - Create Key Pair and Download pem file](./aws-screenshot-5.png)
5. Now Launch Instance. Then you will see you newly created instance 
6. Now you are able to use your pem file to connect to your EC2 instance. [AWS Console - Connect to EC2](./aws-screenshot-6.png)