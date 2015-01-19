################################################################
## Makefile for building Pythia-ROOT-Fastjet code (jetsorter) ## 
## Hannu Siikonen, 12.01.2015                                 ##  
################################################################

SHELL = /bin/bash
CXXFLAGS := -std=c++11 -O2 -ansi -pedantic -W -Wall -Wshadow -Wno-shadow -m64

ROOT=$(shell root-config --cflags --glibs)
PYTHIA=$(shell pythia8-config --cflags --libs)
FASTJET=$(shell fastjet-config --cxxflags --libs)

PYTHIALIB=/usr/local/lib/libpythia8.so

.SECONDEXPANSION:
.PHONY: all clean

all:
	@echo "No make for all objects provided."

parts2tree: $$@.cc $(PYTHIALIB) MinimalEvent.h ../JetSorter/jetsorter_auxiliary.h
	$(ROOTSYS)/bin/rootcint -f $@Dct.cc -c -I. MinimalEvent.h $@LinkDef.h
	$(CXX) $@Dct.cc $@.cc $(PYTHIALIB) -w $(CXXFLAGS) $(ROOT) $(PYTHIA) -ldl -o $@.exe

jetsorttree: $$@.cc
	$(CXX) parts2treeDct.cc $@.cc -w $(CXXFLAGS) $(ROOT) $(FASTJET) -ldl -o $@.exe

modern_jetsorter: $$@.cc
	$(CXX) parts2treeDct.cc $@.cc -w $(CXXFLAGS) $(ROOT) $(FASTJET) -ldl -o $@.exe

testsorter: $$@.cc TreeIterator.o
	$(CXX) parts2treeDct.cc TreeIterator.o $@.cc -w $(CXXFLAGS) $(ROOT) $(FASTJET) -ldl -o $@.exe

TreeIterator.o: TreeIterator.cc
	$(CXX) -c $^ -w $(CXXFLAGS) $(ROOT) $(FASTJET) 


# Clean.
clean:
	rm -f *.log *.exe *.d *.so
