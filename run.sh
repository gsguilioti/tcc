#!/bin/bash

if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <instance_name>"
    exit 1
fi

INSTANCE_NAME=$1

cd ./algorithms || { echo "Error: Could not change to the algorithms directory"; exit 1; }

EXECUTABLES=("2opt" "2optpar2" "3opt" "3optpar2")

for EXEC in "${EXECUTABLES[@]}"; do
    if [ ! -x "$EXEC" ]; then
        echo "Error: Executable $EXEC not found or not executable"
        continue
    fi

    echo "Running $EXEC with instance $INSTANCE_NAME..."
    ./$EXEC "$INSTANCE_NAME"
    echo "Finished running $EXEC"
    echo "----------------------------------------"
done
