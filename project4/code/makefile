# Aim to fix this

CXX = g++

executable = .out # or .exe
detected_OS := $(shell sh -c 'uname 2>/dev/null || echo Unknown')

# Options:
PARALLEL ?= no
OPTIMISE ?= no
ARMA ?= yes
SRC ?= yes 	# are these unnecessary?
EXE ?= yes  # are these unnecessary?

# Flags:
EXTRA_FLAGS =  # dunno
LDFLAGS =  # linker




# Add necessary flag for Mac:
ifeq ($(detected_OS), Darwin)
	CXXFLAGS += -std=c++11
endif

# Use .exe instead of :
ifeq  ($(strip $(EXE)),yes)
    executable = .exe
endif

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

# Add OpenMP flags:
ifeq ($(strip $(PARALLEL)), yes)
	ifeq ($(detected_OS), Darwin)
		EXTRA_FLAGS += -Xpreprocessor
	endif
	EXTRA_FLAGS += -fopenmp
	LDFLAGS += -lomp
endif


all: compile link run

compile:
	$(CXX) $(CXXFLAGS) -c main.cpp $(EXTRA_FLAGS)

link:
	$(CXX) -o main.$(executable) *.o $(LDFLAGS) 

run:
	./main.$(executable) 

debug_valgrind:
	$(CXX) $(CXXFLAGS) -gmain.cpp $(EXTRA_FLAGS) -o  main.$(executable) $(LDFLAGS) && valgrind ./main.$(executable) 

debug_gdb:
	$(CXX) $(CXXFLAGS) -gmain.cpp $(EXTRA_FLAGS) -o  main.$(executable) -$(LDFLAGS) && gdb --args ./main.$(executable)

plots:
	cd plot && python3 analysis.py && cd ..


clean:
	rm -f *.c *.o *.exe *.out