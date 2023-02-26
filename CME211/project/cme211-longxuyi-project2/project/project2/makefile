CXX:= g++
CXXFLAGS:= -std=c++11 -Wall -Wextra -Wconversion
CPPF:= main.cpp CGSolver.cpp COO2CSR.cpp matvecops.cpp heat.cpp sparse.cpp
TARGET:= main
OBJS:= main.o CGSolver.o COO2CSR.o matvecops.o heat.o sparse.o
HPPF:= CGSolver.hpp COO2CSR.hpp matvecops.hpp heat.hpp sparse.hpp

$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) $(OBJS) -o $(TARGET)
%.o: %.cpp $(HPPF)
	$(CXX) -c -o $@ $< $(CXXFLAGS)

clean:
	rm -f main.o $(TARGET) $(OBJS) *-