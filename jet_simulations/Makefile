###################################################################
## Makefile for building ROOT-fastjet analysis of simulated data ## 
## Hannu Siikonen, 04.02.2015                                    ##  
###################################################################

-include config.mk

PYTHIALIB=/usr/local/lib/libpythia8.so

.SECONDEXPANSION:
.PHONY: all clean

INCLUDE := $(INCLUDE) -I.

all:
	@echo "No make for all objects provided."

# Main targets

rootjetsort: src/main_$$@.cpp obj/rootjetsort.o obj/libSimEvent.so
	$(CXX) $^ -w $(INCLUDE) $(CXXFLAGS) $(ROOT) $(FASTJET) -I. -ldl -o $@.exe

ptcut: src/main_$$@.cpp obj/ptcut.o obj/libSimEvent.so
	$(CXX) $^ -w $(INCLUDE) $(CXXFLAGS) $(ROOT) $(FASTJET) -I. -ldl -o $@.exe

obj/SimDct.cc: include/SimEvent.h include/SimLinkDef.h
	$(ROOTSYS)/bin/rootcint -f $@ -c -p $(INCLUDE) $(CFLAGS) $^

# Shared objects

obj/libSimEvent.so: obj/SimEvent.o obj/SimEventDict.o
	$(CXX) -fPIC -shared $(CXXFLAGS) $^ -o $@

# Objects

obj/SimEvent.o: src/SimEvent.cpp
	$(CXX) -fPIC -pthread $^ $(INCLUDE) $(CXXFLAGS) $(ROOT) -c -w -ldl -o $@

obj/SimEventDict.o: obj/SimDct.cc
	$(CXX) $(INCLUDE) $(ROOT) $(CXXFLAGS) -fPIC -pthread -c $^ -o $@

obj/rootjetsort.o: src/RootJetSort.cpp
	$(CXX) -c $^ $(INCLUDE) $(CXXFLAGS) -w $(ROOT) -ldl -o $@

obj/ptcut.o: src/PTCut.cpp
	$(CXX) -c $^ $(INCLUDE) $(CXXFLAGS) -w $(ROOT) -ldl -o $@

# Clean.
clean:
	rm -f *.log *.exe *.d *.so
