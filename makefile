# CC and CFLAGS are varilables
CC = g++
CFLAGS = -c
AR = ar
ARFLAGS = rcv
# -c option ask g++ to compile the source files, but do not link.
# -g option is for debugging version
# -O2 option is for optimized version
OPTFLAGS = -O2

all	: bin/ham_cycle_sat
	@echo -n ""

# optimized version
bin/ham_cycle_sat	: main_opt.o proof.o solver.o file.o
			$(CC) $(OPTFLAGS) main_opt.o proof.o solver.o file.o -o bin/ham_cycle_sat
main_opt.o 	: src/main.cpp
			$(CC) $(CFLAGS) $< -Ilib -o $@
proof.o		: src/Proof.cpp
			$(CC) $(CFLAGS) $(OPTFLAGS) $< -o $@
solver.o	: src/Solver.cpp 
			$(CC) $(CFLAGS) $(OPTFLAGS) $< -o $@
file.o		: src/File.cpp
			$(CC) $(CFLAGS) $(OPTFLAGS) $< -o $@

# clean all the .o and executable files
clean:
		rm -rf *.o lib/*.a lib/*.o bin/*

