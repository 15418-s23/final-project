# Variables
CXX = clang++
CXXFLAGS = -std=c++17 -O3 -Wall -Wextra
SRC_DIR = src
LIB_DIR = $(SRC_DIR)/lib


# Source files
SRCS = $(wildcard $(SRC_DIR)/*.cpp)


# Build rules
all: mcd-sequential mcd-parallel
debug: mcd-sequential-debug mcd-parallel-debug

mcd-sequential:
	$(CXX) $(CXXFLAGS) $(CONFIG) $(SRCS) -o $@

mcd-sequential-debug:
	$(CXX) $(CXXFLAGS) $(CONFIG) -D DEBUG $(SRCS) -o $@

mcd-parallel:
	$(CXX) $(CXXFLAGS) $(CONFIG) -D PARALLEL $(SRCS) -o $@

mcd-parallel-debug:
	$(CXX) $(CXXFLAGS) $(CONFIG) -D DEBUG -D PARALLEL $(SRCS) -o $@


# Clean up
clean:
	rm -rf mcd-sequential mcd-parallel mcd-sequential-debug mcd-parallel-debug


# Phony targets
.PHONY: all debug clean