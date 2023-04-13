EXECUTABLE := mcd # Mesh Collision Detection, not McDonalds...


SOURCES := src/*.cpp
HEADERS := src/*.h


CXX := clang++
CXXFLAGS := -O3 -Wall -g


all: $(EXECUTABLE)

$(EXECUTABLE): 
	$(CXX) -o $@ $(CXXFLAGS) $(SOURCES) $(HEADERS)