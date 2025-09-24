CXX = g++
CXXFLAGS = -std=c++17 -Wall -Wextra -O2 `sdl2-config --cflags`
LDFLAGS = `sdl2-config --libs`


SRC = src/main.cpp


TARGET = bin/game


$(shell mkdir -p bin)

all: $(TARGET)


$(TARGET): $(SRC)
	$(CXX) $(CXXFLAGS) -o $@ $< $(LDFLAGS)

clean:
	rm -f $(TARGET)

.PHONY: all clean
