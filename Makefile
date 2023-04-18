# Variables
CXX = clang++
CXXFLAGS = -std=c++17 -O3 -Wall -Wextra
SRC_DIR = src
LIB_DIR = $(SRC_DIR)/lib


# Source files
SRCS = $(wildcard $(SRC_DIR)/*.cpp)


# Build rules
all: mcd-sequential mcd-parallel

mcd-sequential:
	$(CXX) $(CXXFLAGS) $(SRCS) -o $@

mcd-parallel:
	$(CXX) $(CXXFLAGS) -D PARALLEL $(SRCS) -o $@


# Clean up
clean:
	rm -rf mcd-sequential mcd-parallel


# Phony targets
.PHONY: all clean