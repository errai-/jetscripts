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

pythiatree: src/main_$$@.cpp $(PYTHIALIB) ../obj/minimalevent.so include/jetsorter_auxiliary.h ../obj/PTDct.cc
	$(CXX) $^ $(PYTHIALIB) -w $(INCLUDE) $(CXXFLAGS) $(ROOT) $(PYTHIA) -ldl -o $@.exe

## OBJECTS

../obj/PTDct.cc: ../include/MinimalEvent.h ../include/PTLinkDef.h
	cd .. && $(MAKE) obj/PTDct.cc

../obj/minimalevent.so: ../src/MinimalEvent.cpp
	cd .. && $(MAKE) obj/minimalevent.so

# Clean.
clean:
	rm -f *.log *.exe *.d *.so
