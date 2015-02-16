################################################################
## Makefile for building Pythia-ROOT-Fastjet code (jetsorter) ## 
## Hannu Siikonen, 10.02.2015                                 ##  
################################################################

-include ../config.mk

.SECONDEXPANSION:
.PHONY: all clean

INCLUDE := $(INCLUDE) -I../include -I..

all:
	@echo "No make for all objects provided."

# Main targets

pythia8: src/main_$$@.cpp $(PYTHIALIB) ../obj/libSimEvent.so
	$(CXX) $^ -w $(INCLUDE) $(CXXFLAGS) $(ROOT) $(PYTHIA) -ldl -o $@.exe

test: src/main_$$@.cpp $(PYTHIALIB)
	$(CXX) $^ -w $(INCLUDE) $(CXXFLAGS) $(ROOT) $(PYTHIA) -ldl -o $@.exe

../obj/SimDct.cc: ../include/SimEvent.h ../include/SimLinkDef.h
	cd .. && $(MAKE) obj/SimDct.cc

# Objects

../obj/SimEvent.o: ../src/SimEvent.cpp
	cd .. && $(MAKE) obj/SimEvent.o

../obj/SimEventDict.o: ../obj/SimDct.cc
	cd .. && $(MAKE) obj/SimEventDict.o

# Shared objects

../obj/libSimEvent.so: ../obj/SimEventDict.o ../obj/SimEvent.o
	cd .. && $(MAKE) obj/libSimEvent.so

# Clean.
clean:
	rm -f *.log *.exe *.d *.so

