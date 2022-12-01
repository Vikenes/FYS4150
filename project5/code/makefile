# Makeile for FYS4150 project5

all: compile link run

compile:
	g++ -c -std=c++11 main.cpp src/*.cpp -I include 

link:
	g++ -o main.exe *.o -larmadillo

run: 
	./main.exe 

plots:
	cd plot && python3 analysis.py && cd ..
