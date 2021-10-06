CXX=g++
FLAGS=-Wall -g -O3 -std=c++17
EXE=sexsel_bacteria.exe
CPP=sexsel_bacteria.cpp 


$(EXE) : $(CPP) $(HPP)
	$(CXX) $(FLAGS) -o $(EXE) $(CPP)


.PHONY : clean

clean :
	rm -rf $(EXE)

