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

debug_valgrind:
	g++ -g -std=c++11 main.cpp src/*.cpp -I include -o  main.exe -larmadillo && valgrind ./main.exe 

debug_gdb:
	g++ -g -std=c++11 main.cpp src/*.cpp -I include -o  main.exe -larmadillo && gdb --args ./main.exe

plots:
	cd plot && python3 analysis.py && cd ..
