#!/bin/bash

arg1=$1

gcc -o out motionsim.c -lm

if [ $? -eq 0 ]; then
    echo "Compilation successful. Running the program..."
    ./out "$arg1" 
else
    echo "Compilation failed."
    exit 1 
fi

python3 motionplot.py

