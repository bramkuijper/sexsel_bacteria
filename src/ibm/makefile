CXX=g++
FLAGS=-Wall -O3 -std=c++17 -ggdb
EXE=fishersexsel_bacteria.exe
CPP=fishersexsel_bacteria.cpp individual.cpp
HPP=individual.hpp


$(EXE) : $(CPP) $(HPP)
	$(CXX) $(FLAGS) -o $(EXE) $(CPP)


.PHONY : clean

clean :
	rm -rf $(EXE)

