#!/bin/bash

if [ ! -f "$IDENTITY_FILE" ]; then
    echo "Identity file does not exist"
    exit 1
fi

if [ "$DESTINATION_URL" == "" ]; then
    echo "Destination URL file not set."
    exit 1
fi

if [ "$DESTINATION_PATH" == "" ]; then
    echo "Destination PATH file not set."
    exit 1
fi

echo "start copying to EC2 for all docker image jar files...."

scp -i $IDENTITY_FILE matdb.jar "$DESTINATION_URL:$DESTINATION_PATH"
scp -i $IDENTITY_FILE mlip4tracy_stable.jar "$DESTINATION_URL:$DESTINATION_PATH"
scp -i $IDENTITY_FILE pslibrary4tracy.jar "$DESTINATION_URL:$DESTINATION_PATH"
scp -i $IDENTITY_FILE dft_qe4tracy.jar "$DESTINATION_URL:$DESTINATION_PATH"
scp -i $IDENTITY_FILE ubuntu4tracy.jar "$DESTINATION_URL:$DESTINATION_PATH"
scp -i $IDENTITY_FILE load_docker_images_aws.sh "$DESTINATION_URL:$DESTINATION_PATH"
