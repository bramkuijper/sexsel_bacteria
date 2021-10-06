EXE=solve_fisher.exe
CPP=solve_fisher.cpp main.cpp
HPP=solve_fisher.hpp
CXX=g++
CXXFLAGS=-ggdb -std=c++14 -Wall -O3

all : $(EXE)

$(EXE) : $(CPP) $(HPP)
	$(CXX) $(CXXFLAGS) -o $(EXE) $(CPP)

clean :
	rm -rf $(EXE)
