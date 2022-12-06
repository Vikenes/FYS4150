# Makeile for FYS4150 project5

all: compile link run
all2: compile2 link run
all3: compile3 link run


compile:
	g++ -c -std=c++11 main.cpp src/*.cpp -I include 

compile2:
	g++ -c -std=c++11 main.cpp src/*.cpp -I include -O2

compile3:
	g++ -c -std=c++11 main.cpp src/*.cpp -I include -O3

link:
	g++ -o main.exe *.o -larmadillo

run: 
	./main.exe 

plots:
	cd plot && python3 analysis.py && cd ..




# Clean folder for .c, .o, .exe, .out and .o.tmp files: (should all of these be included?)
clean:
	rm -f *.c *.o *.exe *.out *.o.tmp
