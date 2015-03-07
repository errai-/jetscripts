################################################################
## Makefile for building Pythia-ROOT-Fastjet code (jetsorter) ## 
## Hannu Siikonen, 07.03.2015                                 ##  
################################################################

-include ../config.mk

.SECONDEXPANSION:
.PHONY: all clean

INCLUDE := $(INCLUDE) -I../include -I..

all:
	@echo "No make for all objects provided."

# Main targets

pythia8: src/main_$$@.cpp $(PYTHIALIB) ../obj/libPrtclEvent.so
	$(CXX) $^ -w $(INCLUDE) $(CXXFLAGS) $(ROOT) $(PYTHIA) -ldl -o $@.exe

test: src/main_$$@.cpp $(PYTHIALIB)
	$(CXX) $^ -w $(INCLUDE) $(CXXFLAGS) $(ROOT) $(PYTHIA) -ldl -o $@.exe

../obj/PrtclDct.cc: ../include/PrtclEvent.h ../include/PrtclLinkDef.h
	cd .. && $(MAKE) obj/PrtclDct.cc

# Objects

../obj/PrtclEvent.o: ../src/PrtclEvent.cpp
	cd .. && $(MAKE) obj/PrtclEvent.o

../obj/PrtclEventDict.o: ../obj/PrtclDct.cc
	cd .. && $(MAKE) obj/PrtclEventDict.o

# Shared objects

../obj/libPrtclEvent.so: ../obj/PrtclEventDict.o ../obj/PrtclEvent.o
	cd .. && $(MAKE) obj/libPrtclEvent.so

# Clean.
clean:
	rm -f *.log *.exe *.d *.so

