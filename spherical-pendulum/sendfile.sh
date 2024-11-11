#!/bin/bash

if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <file-to-copy>"
    exit 1
fi

LOCAL_FILE=$1
REMOTE_USER="richardanderson"
REMOTE_HOST="192.168.1.177"
REMOTE_PATH="/Users/richardanderson/remote-work/"
PASSWORD="annaevin1214"

sshpass -p "$PASSWORD" scp "$LOCAL_FILE" "$REMOTE_USER@$REMOTE_HOST:$REMOTE_PATH"

if [ $? -eq 0 ]; then
    echo "File successfully copied to $REMOTE_USER@$REMOTE_HOST:$REMOTE_PATH"
else
    echo "Error occurred while copying the file."
fi

