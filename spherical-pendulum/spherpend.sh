#!/bin/bash

gcc -o pendout spherpendulum.c

if [ $? -eq 0 ]; then
    echo "Compilation successful. Running the program..."
    ./pendout
else
    echo "Compilation failed."
fi

python3 sphervisp.py

