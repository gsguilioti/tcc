#!/bin/bash

if [ "$#" -lt 1 ]; then
    echo "Usage: $0 <instance1> [instance2] ... [instanceN]"
    exit 1
fi

cd ./algorithms || { echo "Error: Could not change to the algorithms directory"; exit 1; }

EXECUTABLES=("2opt" "2optpar2" "3opt" "3optpar2")

for INSTANCE_NAME in "$@"; do
    echo "Processing instance: $INSTANCE_NAME"
    
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
done
