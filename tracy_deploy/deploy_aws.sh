IDENTITY_FILE="matdb-test.pem" 
if [ -f "$IDENTITY_FILE" ]; then
    echo "start copying to EC2 for all docker image jar files...."
else 
    echo ".pem file does not exist"
    exit 1
fi
scp -i $IDENTITY_FILE  matdb_stable.jar ubuntu@ec2-18-220-253-111.us-east-2.compute.amazonaws.com:~/
scp -i $IDENTITY_FILE  mlip4tracy_stable.jar ubuntu@ec2-18-220-253-111.us-east-2.compute.amazonaws.com:~/
scp -i $IDENTITY_FILE  tracy_science.jar ubuntu@ec2-18-220-253-111.us-east-2.compute.amazonaws.com:~/
scp -i $IDENTITY_FILE  pslibrary4tracy.jar ubuntu@ec2-18-220-253-111.us-east-2.compute.amazonaws.com:~/
scp -i $IDENTITY_FILE  dft_qe4tracy.jar ubuntu@ec2-18-220-253-111.us-east-2.compute.amazonaws.com:~/
scp -i $IDENTITY_FILE  ubuntu4tracy.jar ubuntu@ec2-18-220-253-111.us-east-2.compute.amazonaws.com:~/
scp -i $IDENTITY_FILE  ubuntu.jar ubuntu@ec2-18-220-253-111.us-east-2.compute.amazonaws.com:~/
scp -i $IDENTITY_FILE load_docker_images_aws.sh ubuntu@ec2-18-220-253-111.us-east-2.compute.amazonaws.com:~/
