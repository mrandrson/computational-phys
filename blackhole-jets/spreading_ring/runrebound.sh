#!/bin/bash

echo "Script to run:"
scriptname=problem

gcc "$scriptname.c" -o "$scriptname" -I/Users/richardanderson/rebound/rebound -I/opt/homebrew/opt/hdf5/include -L/opt/homebrew/opt/hdf5/lib -lhdf5 ./librebound.so -lm

if [ $? -eq 0 ]; then
    ./"$scriptname"
else
    echo "Compilation failed."
fi

