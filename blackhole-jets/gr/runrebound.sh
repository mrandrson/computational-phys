#!/bin/bash

echo "Script to run:"
scriptname=problem

export DYLD_LIBRARY_PATH="/Users/richardanderson/reboundx:$DYLD_LIBRARY_PATH"

gcc "$scriptname.c" -o "$scriptname" \
    -I/Users/richardanderson/rebound/rebound \
    -I/Users/richardanderson/reboundx/reboundx \
    -I/opt/homebrew/opt/hdf5/include \
    -L/opt/homebrew/opt/hdf5/lib -lhdf5 \
    -L/Users/richardanderson/reboundx -lreboundx \
    ./librebound.so -lm

if [ $? -eq 0 ]; then
    ./"$scriptname"
else
    echo "Compilation failed."
fi

