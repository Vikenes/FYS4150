# Makeile for FYS4150 project5


### GENERAL STUFF

CXX = g++
CXXFLAGS += -std=c++11

executable = .exe # or .out
detected_OS := $(shell sh -c 'uname 2>/dev/null || echo Unknown') # will be Darwin for Mac OS X

# Options:
PARALLEL ?= no
OPTIMISE ?= no
ARMA ?= yes
SRC ?= yes 	# are these unnecessary?

# Flags:
EXTRA_FLAGS =  # dunno
LDFLAGS =  # linker

# # Add necessary flag for Mac:
# ifeq ($(detected_OS), Darwin)
# 	CXXFLAGS += -std=c++11
# endif

# Add larmadillo flag:
ifeq  ($(strip $(ARMA)),yes)
    LDFLAGS += -larmadillo
endif

# For header files:
ifeq  ($(strip $(SRC)),yes)
    EXTRA_FLAGS += src/*.cpp -I include
endif

# Add O3 flag:
ifeq  ($(strip $(OPTIMISE)),yes)
    CXXFLAGS += -O3
endif

# # Add OpenMP flags:
# ifeq ($(strip $(PARALLEL)), yes)
# 	ifeq ($(detected_OS), Darwin)
# 		EXTRA_FLAGS += -Xpreprocessor
# 	endif
# 	EXTRA_FLAGS += -fopenmp
# 	LDFLAGS += -lomp
# endif

WHICH ?= ALL # arg. to main.exe


# General run command: 
all: compile link run 

# Compile main.cpp:
compile:
	$(CXX) $(CXXFLAGS) -c main.cpp $(EXTRA_FLAGS)

# Link main.exe:
link:
	$(CXX) -o main$(executable) *.o $(LDFLAGS) 

# Run main.exe: + project specific arg.
run:
	./main$(executable) $(WHICH)

# Debug, not sure:
debug_valgrind:
	$(CXX) $(CXXFLAGS) -gmain.cpp $(EXTRA_FLAGS) -o  main$(executable) $(LDFLAGS) && valgrind ./main$(executable) 

# Debuh, not sure:
debug_gdb:
	$(CXX) $(CXXFLAGS) -gmain.cpp $(EXTRA_FLAGS) -o  main$(executable) -$(LDFLAGS) && gdb --args ./main$(executable)



# Create and show all plots:
plots:
	cd plot && python3 analysis.py $(WHICH) figures && cd ..
# Aiming to have options here!


# Create and show all animations:
animations:
	cd plot && python3 analysis.py $(WHICH) animate && cd ..



# Clean folder for .c, .o, .exe, .out and .o.tmp files: (should all of these be included?)
clean:
	rm -f *.c *.o *.exe *.out *.o.tmp



### OLD


# all: compile link run
# all2: compile2 link run
# all3: compile3 link run


# compile:
# 	g++ -c -std=c++11 main.cpp src/*.cpp -I include 

# compile2:
# 	g++ -c -std=c++11 main.cpp src/*.cpp -I include -O2

# compile3:
# 	g++ -c -std=c++11 main.cpp src/*.cpp -I include -O3

# link:
# 	g++ -o main.exe *.o -larmadillo

# run: 
# 	./main.exe 

# plots:
# 	cd plot && python3 analysis.py && cd ..




# # Clean folder for .c, .o, .exe, .out and .o.tmp files: (should all of these be included?)
# clean:
# 	rm -f *.c *.o *.exe *.out *.o.tmp
