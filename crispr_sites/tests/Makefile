PROGRAM_NAME=crispr_sites
PROGRAM_VERSION := $(shell git describe --dirty --always --tags)

CPPFLAGS=--std=c++11 -O3 

tests_all : main.o scan_stdin.o
	g++ $(CPPFLAGS) -o tests_all main.o scan_stdin.o crispr_sites.o

main.o scan_stdin.o crispr_sites.o : main.cpp scan_stdin.cpp ../crispr_sites.cpp
	g++ $(CPPFLAGS) -DPROGRAM_VERSION=\"$(PROGRAM_VERSION)\" -DPROGRAM_NAME=\"$(PROGRAM_NAME)\" -c main.cpp scan_stdin.cpp -DUNIT_TESTS ../crispr_sites.cpp

.PHONY: clean

clean:
	rm -f tests_all *.o

