################################################################
## Makefile for building Pythia-ROOT-Fastjet code (jetsorter) ## 
## Hannu Siikonen, 12.01.2015                                 ##  
################################################################

SHELL = /bin/bash
CXXFLAGS := -O2 -ansi -pedantic -W -Wall -Wshadow -Wno-shadow -m64 -std=c++11  

ROOT=$(shell root-config --cflags --glibs)
PYTHIA=$(shell pythia8-config --cflags --libs)
FASTJET=$(shell fastjet-config --cxxflags --libs)

PYTHIALIB=/usr/local/lib/libpythia8.so

.SECONDEXPANSION:
.PHONY: all clean

all:
	@echo "No make for all objects provided."

# Main targets

jetsorter_modern ptcut_naive_modern: src/$$@.cc $(PYTHIALIB) obj/minimalevent.o
	$(CXX) $^ -w $(CXXFLAGS) $(ROOT) $(PYTHIA) $(FASTJET) -ldl -o $@.exe 

jetsorter: src/main_$$@.cpp obj/jetsorter.o obj/minimalevent.o
	$(CXX) $^ $(CXXFLAGS) -w $(ROOT) $(FASTJET) -ldl -o $@.exe

pythiasaver: src/main_$$@.cpp obj/pythiasaver.o obj/minimalevent.o
	$(CXX) $^ $(CXXFLAGS) -w $(ROOT) $(PYTHIA) -ldl -o $@.exe

## OBJECTS
	
obj/jetsorter.o: src/JetSorter.cpp $(PYTHIALIB)
	$(CXX) -c $^ $(CXXFLAGS) -w $(ROOT) $(FASTJET) -ldl -o $@
	
obj/pythiasaver.o: src/PythiaSaver.cpp $(PYTHIALIB)
	$(CXX) -c $^ $(CXXFLAGS) -w $(ROOT) $(PYTHIA) -ldl -o $@

obj/minimalevent.o: src/MinimalEvent.cpp
	$(CXX) -c $^ $(CXXFLAGS) -w -ldl -o obj/minimalevent.o
	
# Clean.
clean:
	rm -f *.log *.exe *.d *.so
