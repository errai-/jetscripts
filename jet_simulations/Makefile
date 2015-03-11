###################################################################
## Makefile for building ROOT-fastjet analysis of simulated data ## 
## Hannu Siikonen, 11.03.2015                                    ##  
###################################################################

-include config.mk

PYTHIALIB=/usr/local/lib/libpythia8.so

.SECONDEXPANSION:
.PHONY: all clean

INCLUDE := $(INCLUDE) -I.

all:
	@echo "No make for all objects provided."

# Main targets

main_jetanalysis: jetanalysis/$$@.cpp obj/jetanalysis.o lib/libPrtclEvent.so lib/libJetEvent.so
	$(CXX) $^ -w $(INCLUDE) $(CXXFLAGS) $(ROOT) $(FASTJET) -I. -ldl -o \
jetanalysis.exe

rootjetsort: src/main_$$@.cpp obj/rootjetsort.o obj/libPrtclEvent.so
	$(CXX) $^ -w $(INCLUDE) $(CXXFLAGS) $(ROOT) $(FASTJET) -I. -ldl -o $@.exe

ptcut: src/main_$$@.cpp obj/ptcut.o obj/libPrtclEvent.so
	$(CXX) $^ -w $(INCLUDE) $(CXXFLAGS) $(ROOT) $(FASTJET) -I. -ldl -o $@.exe

obj/PrtclDct.cc: events/PrtclEvent.h events/PrtclLinkDef.h
	$(ROOTSYS)/bin/rootcint -f $@ -c -p $(INCLUDE) $(CFLAGS) $^

obj/JetDct.cc: events/JetEvent.h events/JetLinkDef.h
	$(ROOTSYS)/bin/rootcint -f $@ -c -p $(INCLUDE) $(CFLAGS) $^

main_pythia8: pythia8/$$@.cpp $(PYTHIALIB) lib/libPrtclEvent.so
	$(CXX) $^ -w $(INCLUDE) $(CXXFLAGS) $(ROOT) $(PYTHIA) -ldl -o \
pythia8.exe

main_test: pythia8/$$@.cpp $(PYTHIALIB)
	$(CXX) $^ -w $(INCLUDE) $(CXXFLAGS) $(ROOT) $(PYTHIA) -ldl -o test.exe

# Shared objects

lib/libPrtclEvent.so: obj/PrtclEvent.o obj/PrtclEventDict.o
	$(CXX) -fPIC -shared $(CXXFLAGS) $^ -o $@

lib/libJetEvent.so: obj/JetEvent.o obj/JetEventDict.o
	$(CXX) -fPIC -shared $(CXXFLAGS) $^ -o $@

# Objects

obj/PrtclEvent.o: events/PrtclEvent.cpp
	$(CXX) -fPIC -pthread $^ $(INCLUDE) $(CXXFLAGS) $(ROOT) -c -w -ldl -o $@

obj/PrtclEventDict.o: obj/PrtclDct.cc
	$(CXX) $(INCLUDE) $(ROOT) $(CXXFLAGS) -fPIC -pthread -c $^ -o $@

obj/JetEvent.o: events/JetEvent.cpp
	$(CXX) -fPIC -pthread $^ $(INCLUDE) $(CXXFLAGS) $(ROOT) -c -w -ldl -o $@

obj/JetEventDict.o: obj/JetDct.cc
	$(CXX) $(INCLUDE) $(ROOT) $(CXXFLAGS) -fPIC -pthread -c $^ -o $@

obj/rootjetsort.o: src/RootJetSort.cpp
	$(CXX) -c $^ $(INCLUDE) $(CXXFLAGS) -w $(ROOT) -ldl -o $@

obj/jetanalysis.o: jetanalysis/JetAnalysis.cpp lib/libJetEvent.so lib/libPrtclEvent.so
	$(CXX) -c $^ $(INCLUDE) $(CXXFLAGS) -w $(ROOT) -ldl -o $@

obj/ptcut.o: src/PTCut.cpp
	$(CXX) -c $^ $(INCLUDE) $(CXXFLAGS) -w $(ROOT) -ldl -o $@

# Clean.
clean:
	rm -f *.log *.exe *.d *.so
