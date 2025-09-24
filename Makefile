# Compiler
CXX = g++
CXXFLAGS = -std=c++17 -Wall -Wextra -O2 `sdl2-config --cflags`
LDFLAGS = `sdl2-config --libs`

# Source file
SRC = src/main.cpp

# Output binary
TARGET = bin/game

# Create bin directory if it doesn't exist1
$(shell mkdir -p bin)

# Default target
all: $(TARGET)

# Compile
$(TARGET): $(SRC)
	$(CXX) $(CXXFLAGS) -o $@ $< $(LDFLAGS)

# Clean
clean:
	rm -f $(TARGET)

.PHONY: all clean
