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

jetanalysis: src/main_$$@.cpp obj/jetanalysis.o obj/libPrtclEvent.so obj/libJetEvent.so
	$(CXX) $^ -w $(INCLUDE) $(CXXFLAGS) $(ROOT) $(FASTJET) -I. -ldl -o $@.exe

rootjetsort: src/main_$$@.cpp obj/rootjetsort.o obj/libPrtclEvent.so
	$(CXX) $^ -w $(INCLUDE) $(CXXFLAGS) $(ROOT) $(FASTJET) -I. -ldl -o $@.exe

ptcut: src/main_$$@.cpp obj/ptcut.o obj/libPrtclEvent.so
	$(CXX) $^ -w $(INCLUDE) $(CXXFLAGS) $(ROOT) $(FASTJET) -I. -ldl -o $@.exe

obj/PrtclDct.cc: include/PrtclEvent.h include/PrtclLinkDef.h
	$(ROOTSYS)/bin/rootcint -f $@ -c -p $(INCLUDE) $(CFLAGS) $^

obj/JetDct.cc: include/JetEvent.h include/JetLinkDef.h
	$(ROOTSYS)/bin/rootcint -f $@ -c -p $(INCLUDE) $(CFLAGS) $^

# Shared objects

obj/libPrtclEvent.so: obj/PrtclEvent.o obj/PrtclEventDict.o
	$(CXX) -fPIC -shared $(CXXFLAGS) $^ -o $@

obj/libJetEvent.so: obj/JetEvent.o obj/JetEventDict.o
	$(CXX) -fPIC -shared $(CXXFLAGS) $^ -o $@

# Objects

obj/PrtclEvent.o: src/PrtclEvent.cpp
	$(CXX) -fPIC -pthread $^ $(INCLUDE) $(CXXFLAGS) $(ROOT) -c -w -ldl -o $@

obj/PrtclEventDict.o: obj/PrtclDct.cc
	$(CXX) $(INCLUDE) $(ROOT) $(CXXFLAGS) -fPIC -pthread -c $^ -o $@

obj/JetEvent.o: src/JetEvent.cpp
	$(CXX) -fPIC -pthread $^ $(INCLUDE) $(CXXFLAGS) $(ROOT) -c -w -ldl -o $@

obj/JetEventDict.o: obj/JetDct.cc
	$(CXX) $(INCLUDE) $(ROOT) $(CXXFLAGS) -fPIC -pthread -c $^ -o $@

obj/rootjetsort.o: src/RootJetSort.cpp
	$(CXX) -c $^ $(INCLUDE) $(CXXFLAGS) -w $(ROOT) -ldl -o $@

obj/jetanalysis.o: src/JetAnalysis.cpp obj/libJetEvent.so obj/libPrtclEvent.so
	$(CXX) -c $^ $(INCLUDE) $(CXXFLAGS) -w $(ROOT) -ldl -o $@

obj/ptcut.o: src/PTCut.cpp
	$(CXX) -c $^ $(INCLUDE) $(CXXFLAGS) -w $(ROOT) -ldl -o $@

# Clean.
clean:
	rm -f *.log *.exe *.d *.so
