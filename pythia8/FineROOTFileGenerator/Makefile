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

event2tree: $$@.cc $(PYTHIALIB) 
	$(CXX) $^ -w $(CXXFLAGS) $(ROOT) $(PYTHIA) $(FASTJET) -ldl -o $@.exe 

clean:
	rm -f *.log *.exe *.d *.so 
