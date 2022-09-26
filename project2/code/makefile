all: compile link run

compile:
	g++ -c main.cpp src/*.cpp -I include

link:
	g++ -o main.exe *.o -larmadillo

run:
	./main.exe 

