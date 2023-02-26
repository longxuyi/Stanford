CXX:= g++
CXXFLAGS:= -O3 -std=c++11 -Wall -Wextra -Wconversion -Wpedantic
CPPF:= main.cpp hw6.cpp image.cpp
TARGET:= main
OBJS:= main.o hw6.o image.o 
HPPF:= hw6.hpp image.hpp

$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) $(OBJS) -ljpeg -o $(TARGET)
%.o: %.cpp $(HPPF)
	$(CXX) -c -o $@ $< $(CXXFLAGS)

clean:
	rm -f main.o $(TARGET) $(OBJS) *-