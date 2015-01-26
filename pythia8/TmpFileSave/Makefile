################################################################
## Makefile for building Pythia-ROOT-Fastjet code (jetsorter) ## 
## Hannu Siikonen, 12.01.2015                                 ##  
################################################################

SHELL = /bin/bash
CFLAGS := -O2 -ansi -pedantic -W -Wall -Wshadow -Wno-shadow -m64
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

pythiatree: src/main_$$@.cpp $(PYTHIALIB) obj/minimalevent.o ../JetSorter/jetsorter_auxiliary.h obj/PTDct.cc
	$(CXX) $^ $(PYTHIALIB) -w $(CXXFLAGS) $(ROOT) $(PYTHIA) -I. -ldl -o $@.exe

rootjetsort: src/main_$$@.cpp obj/rootjetsort.o obj/jetstore.o obj/minimalevent.o obj/PTDct.cc
	$(CXX) $^ -w $(CXXFLAGS) $(ROOT) $(FASTJET) -I. -ldl -o $@.exe

jetsorttree: src/$$@.cc obj/PTDct.cc obj/minimalevent.o obj/jetstore.o
	$(CXX) $^ -w $(CXXFLAGS) $(ROOT) $(FASTJET) -ldl -I. -o $@.exe

jetsorter_modern ptcut_naive_modern: src/$$@.cc $(PYTHIALIB) obj/minimalevent.o
	$(CXX) $^ -w $(CXXFLAGS) $(ROOT) $(PYTHIA) $(FASTJET) -ldl -o $@.exe 

jetsorter: src/main_$$@.cpp obj/libP2TLib.so obj/minimalevent.o
	$(CXX) $^ $(CXXFLAGS) -w $(ROOT) $(FASTJET) -ldl -o $@.exe

pythiasaver: src/main_$$@.cpp obj/pythiasaver.o obj/minimalevent.o
	$(CXX) $^ $(CXXFLAGS) -w $(ROOT) $(PYTHIA) -ldl -o $@.exe

## OBJECTS

obj/libP2TLib.so: obj/PTDct.cc src/MinimalEvent.cpp
	g++ -shared -o $@ $(ROOT) $(CXXFLAGS) -I. $^ -fPIC

obj/PTDct.cc: include/MinimalEvent.h include/PTLinkDef.h
	$(ROOTSYS)/bin/rootcint -f $@ -c -p $(CFLAGS) -I. $^

obj/rootjetsort.o: src/RootJetSort.cpp
	$(CXX) -c $^ $(CXXFLAGS) -w $(ROOT) -ldl -o $@

obj/jetstore.o: src/JetStore.cpp 
	$(CXX) -c $^ $(CXXFLAGS) -w $(ROOT) $(FASTJET) -ldl -o $@

obj/jetsorter.o: src/JetSorter.cpp $(PYTHIALIB)
	$(CXX) -c $^ $(CXXFLAGS) -w $(ROOT) $(FASTJET) -ldl -o $@
	
obj/pythiasaver.o: src/PythiaSaver.cpp $(PYTHIALIB)
	$(CXX) -c $^ $(CXXFLAGS) -w $(ROOT) $(PYTHIA) -ldl -o $@

obj/minimalevent.o: src/MinimalEvent.cpp
	$(CXX) -c $^ $(CXXFLAGS) -w -ldl -o obj/minimalevent.o
	
# Clean.
clean:
	rm -f *.log *.exe *.d *.so
