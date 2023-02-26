# main: main.cpp CGSolver.cpp COO2CSR.cpp matvecops.cpp CGSolver.hpp COO2CSR.hpp matvecops.hpp
# 	g++ -Wall -Wextra -Wconversion -o main main.cpp CGSolver.cpp COO2CSR.cpp matvecops.cpp
# clean:
# 	rm *.o main

CXX:= g++
CXXFLAGS:= -std=c++11 -Wall -Wextra -Wconversion
CPPF:= main.cpp CGSolver.cpp COO2CSR.cpp matvecops.cpp 
TARGET:= main
OBJS:= main.o CGSolver.o COO2CSR.o matvecops.o 
HPPF:= CGSolver.hpp COO2CSR.hpp matvecops.hpp

$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) $(OBJS) -o $(TARGET)
%.o: %.cpp $(HPPF)
	$(CXX) -c -o $@ $< $(CXXFLAGS)

clean:
	rm -f main.o $(TARGET) $(OBJS) *-