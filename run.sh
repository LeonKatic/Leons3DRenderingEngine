#!/bin/bash

make

if [ $? -eq 0 ]; then
    echo "Compilation successful. Running the game..."
    ./bin/game
else
    echo "Compilation failed."
fi
