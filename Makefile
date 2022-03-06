all: smallpt.cpp
	g++ -O3 -fopenmp smallpt.cpp -o smallpt