CXX=g++
FLAGS=-Wall -g -O3 -std=c++17
EXE=fishersexsel_Gandon.exe
CPP=fishersexsel_bacteria_Gandon_new.cpp individual.cpp
HPP=individual.hpp


$(EXE) : $(CPP) $(HPP)
	$(CXX) $(FLAGS) -o $(EXE) $(CPP)


.PHONY : clean

clean :
	rm -rf $(EXE)

