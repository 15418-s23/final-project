# Variables
TARGET = mcd # Mesh Collision Detection, not McDonalds...
CXX = clang++
CXXFLAGS = -std=c++17 -O3 -Wall -Wextra
SRC_DIR = src
LIB_DIR = $(SRC_DIR)/lib


# Source files
SRCS = $(wildcard $(SRC_DIR)/*.cpp)


# Build rules
all: $(TARGET)

$(TARGET):
	$(CXX) $(CXXFLAGS) $(SRCS) -o $@


# Clean up
clean:
	rm -rf $(TARGET)


# Phony targets
.PHONY: all clean