#############################################################
## Makefile for building Pythia-ROOT-Fastjet files (ptcut) ##
## Hannu Siikonen, 12.01.2014                              ##
#############################################################

# The Makefile might require some fine tuning depending on the
# install location of fastjet and pythia8

SHELL = /bin/bash
LIBRARY_PATH := -Wl,-rpath /usr/local/lib # pythia8 library location
CXXFLAGS := -O2 -ansi -pedantic -W -Wall -Wshadow -Wno-shadow -m64

ROOT=$(shell root-config --cflags --glibs)
PYTHIA=$(shell pythia8-config --cflags --libs)
FASTJET=$(shell fastjet-config --cxxflags --libs)

PYTHIALIB=/usr/local/lib/libpythia8.so

.SECONDEXPANSION:
.PHONY: all clean

all:
	@echo "No make for all objects provided."

ptcut ptcut_progressive ptcut_stable ptcut_stable_sum ptcut_datahists: $$@.cc $(PYTHIALIB) 
	$(CXX) $^ -w $(CXXFLAGS) $(ROOT) $(PYTHIA) $(FASTJET) -ldl -o $@.exe 

clean:
	rm -f $(EXE) *.log *.exe *.d *.so 

