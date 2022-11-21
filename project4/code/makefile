# Aim to fix this

CXX = g++

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



# Add necessary flag for Mac:
ifeq ($(detected_OS), Darwin)
	CXXFLAGS += -std=c++11
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


compile_anal:
	g++ -std=c++11 analytical.cpp src/*.cpp -I include -o analytical.exe -larmadillo  
run_anal:
	./analytical.exe $N $(N_save) $T

compile_equi:
	g++ -std=c++11 equilibriation_time.cpp src/*.cpp -I include -o equilibriation_time.exe -larmadillo 
run_equi:
	./equilibriation_time.exe $N $T $(order)

compile_pdf:
	g++ -std=c++11 prob_distr.cpp src/*.cpp -I include -o prob_distr.exe -larmadillo
run_pdf:
	./prob_distr.exe $(N_samples) $(N_eq_steps) $T $(order)

parallel_compile:
	g++ -std=c++11 -O2 parallel_system.cpp src/*.cpp -I include -o parallel_system.exe -larmadillo -fopenmp
unparallel_compile:
	g++ -std=c++11 -O2 parallel_system.cpp src/*.cpp -I include -o unparallel_system.exe -larmadillo
parallel_run:
	./parallel_system.exe $L $(T0) $(T1) $(nT) $(NMC) $(Neq)
unparallel_run:
	./unparallel_system.exe $L $(T0) $(T1) $(nT) $(NMC) $(Neq)

# General run command:
all: compile link run
# Aiming to have options here!

# Compile main.cpp:
compile:
	$(CXX) $(CXXFLAGS) -c main.cpp $(EXTRA_FLAGS)

# Link main.exe:
link:
	$(CXX) -o main$(executable) *.o $(LDFLAGS) 

# Run main.exe:
run:
	./main$(executable)

# Debug, not sure:
debug_valgrind:
	$(CXX) $(CXXFLAGS) -gmain.cpp $(EXTRA_FLAGS) -o  main$(executable) $(LDFLAGS) && valgrind ./main$(executable) 

# Debuh, not sure:
debug_gdb:
	$(CXX) $(CXXFLAGS) -gmain.cpp $(EXTRA_FLAGS) -o  main$(executable) -$(LDFLAGS) && gdb --args ./main$(executable)

# Create and show all plots:
plots:
	cd plot && python3 analysis.py && cd ..
# Aiming to have options here!

# Clean folder for .c, .o, .exe, .out and .o.tmp files: (should all of these be included?)
clean:
	rm -f *.c *.o *.exe *.out *.o.tmp
