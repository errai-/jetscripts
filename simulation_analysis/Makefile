###################################################################
## Makefile for building ROOT-fastjet analysis of simulated data ## 
## Hannu Siikonen, 04.02.2015                                    ##  
###################################################################

-include config.mk

PYTHIALIB=/usr/local/lib/libpythia8.so

.SECONDEXPANSION:
.PHONY: all clean

all:
	@echo "No make for all objects provided."

# Main targets

rootjetsort: src/main_$$@.cpp obj/rootjetsort.o obj/minimalevent.so obj/PTDct.cc
	$(CXX) $^ -w $(CXXFLAGS) $(ROOT) $(FASTJET) -I. -ldl -o $@.exe

## OBJECTS

obj/minimalevent.so: src/MinimalEvent.cpp
	$(CXX) -shared -fPIC $^ $(INCLUDE) $(CXXFLAGS) $(ROOT) -w -ldl -o obj/minimalevent.so

obj/libP2TLib.so: obj/PTDct.cc src/MinimalEvent.cpp
	g++ -shared -o $@ $(INCLUDE) $(ROOT) $(CXXFLAGS) -I. $^ -fPIC

obj/PTDct.cc: include/MinimalEvent.h include/PTLinkDef.h
	$(ROOTSYS)/bin/rootcint -f $@ -c -p $(INCLUDE) $(CFLAGS) $^

obj/rootjetsort.o: src/RootJetSort.cpp
	$(CXX) -c $^ $(INCLUDE) $(CXXFLAGS) -w $(ROOT) -ldl -o $@

# Clean.
clean:
	rm -f *.log *.exe *.d *.so
