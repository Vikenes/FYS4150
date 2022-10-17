all: compile link run

compile:
	g++ -c -std=c++11 main.cpp src/*.cpp -I include

link:
	g++ -o main.exe *.o -larmadillo

run:
	./main.exe 

debug_valgrind:
	g++ -g -std=c++11 main.cpp src/*.cpp -I include -o  main.exe -larmadillo && valgrind ./main.exe 

debug_gdb:
	g++ -g -std=c++11 main.cpp src/*.cpp -I include -o  main.exe -larmadillo && gdb --args ./main.exe

plots:
	cd plot && python3 plot.py && cd ..
