#!/bin/bash

REMOTE_USER="muti9295"
REMOTE_SERVER="terra.colorado.edu"
REMOTE_DIR="/data/becker/muti9295/PulseReconstruction/runs"
LOCAL_DIR="results/"

FILES=$(ssh "$REMOTE_USER@$REMOTE_SERVER" "find $REMOTE_DIR -type f -name '*.mat' ! -name '*checkpoint*'")

if [ -n "$FILES" ]; then
    for FILE in $FILES; do
         rsync -av --ignore-existing --progress "$REMOTE_USER@$REMOTE_SERVER:$FILE" "$LOCAL_DIR"
    done
else
    echo "No files found."
fi
