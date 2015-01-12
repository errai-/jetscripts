################################################################
## Makefile for building Pythia-ROOT-Fastjet code (jetsorter) ## 
## Hannu Siikonen, 12.01.2015                                 ##  
################################################################

SHELL = /bin/bash
CXXFLAGS := -O2 -ansi -pedantic -W -Wall -Wshadow -Wno-shadow -m64

ROOT=$(shell root-config --cflags --glibs)
PYTHIA=$(shell pythia8-config --cflags --libs)
FASTJET=$(shell fastjet-config --cxxflags --libs)

PYTHIALIB=/usr/local/lib/libpythia8.so

.SECONDEXPANSION:
.PHONY: all clean

all:
	@echo "No make for all objects provided."

# Basic jetsorter
jetsorter: $$@.cc $(PYTHIALIB)
	$(CXX) $^ -w $(CXXFLAGS) $(ROOT) $(PYTHIA) $(FASTJET) -ldl -o $@.exe 

# The main-files are from the pythia8 source, used to test that everything works
main71 main72: $$@.cc $(PYTHIALIB)
	$(CXX) $^ $(CXXFLAGS) $(PYTHIA) $(FASTJET) -o $@.exe 

main91: $$@.cc $(PYTHIALIB) 
	$(CXX) $^ -w $(CXXFLAGS) $(ROOT) $(PYTHIA) -ldl -o $@.exe 

main92: $$@.cc $$@.h $$@LinkDef.h $(PYTHIALIB) 
	`rootcint -f $@Dct.cc -c `pythia8-config --cflags` $@.h $@LinkDef.h`
	$(CXX) $@Dct.cc $^ -w $(CXXFLAGS) $(ROOT) $(PYTHIA) -o $@.exe 

# Clean.
clean:
	rm -f *.log *.exe *.d *.so
