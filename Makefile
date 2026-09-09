CXX ?= g++
CXXFLAGS ?= -O2 -std=c++17 -Wall -Wextra

run: main.cpp
	$(CXX) $(CXXFLAGS) -o run main.cpp

clean:
	rm -f run

.PHONY: clean
