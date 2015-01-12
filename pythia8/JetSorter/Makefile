################################################################
## Makefile for building Pythia-ROOT-Fastjet code (jetsorter) ## 
## Hannu Siikonen, 12.01.2015                                 ##  
################################################################

SHELL = /bin/bash
LIBRARY_PATH := -Wl,-rpath /usr/local/lib # pythia8 library location
CXXFLAGS := -O2 -ansi -pedantic -W -Wall -Wshadow -Wno-shadow -m64

ROOT=$(shell root-config --cflags --glibs)
PYTHIA=$(shell pythia8-config --cflags --libs)
FASTJET=$(shell fastjet-config --cxxflags --libs)
PYTHIALIB=/usr/local/lib/libpythia8.so
################################################################################
# RULES: Definition of the rules used to build the PYTHIA examples.
################################################################################

# Rules without physical targets (secondary expansion for specific rules).
.SECONDEXPANSION:
.PHONY: all clean

# All targets (no default behavior).
all:
	@echo "No make for all objects provided."

# Basic jetsorter
jetsorter: $$@.cc $(PYTHIALIB)
	$(CXX) $@.cc -w $(CXXFLAGS) $(ROOT) $(PYTHIA) $(FASTJET) -ldl -o $@.exe 

main71 main72: $$@.cc $(PYTHIALIB)
	$(CXX) $^ $(CXX_FLAGS) $(PYTHIA) $(FASTJET) -o $@.exe 

main91: $$@.cc $(PYTHIALIB) 
	$(CXX) $^ -w $(CXX_FLAGS) $(ROOT) $(PYTHIA) -ldl -o $@.exe 

main92: $$@.cc $$@.h $$@LinkDef.h $(PYTHIALIB) 
	`rootcint -f $@Dct.cc -c `pythia8-config --cflags` $@.h $@LinkDef.h`
	$(CXX) $@Dct.cc $^ -w $(CXXFLAGS) $(ROOT) $(PYTHIA) -o $@.exe 

# Clean.
clean:
	rm -f *.log *.exe *.d *.so
