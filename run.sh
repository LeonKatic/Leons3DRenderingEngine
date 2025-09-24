#!/bin/bash

# Compile the project
make

# Check if compilation succeeded
if [ $? -eq 0 ]; then
    echo "Compilation successful. Running the game..."
    ./bin/game
else
    echo "Compilation failed."
fi
