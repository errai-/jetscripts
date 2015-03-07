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

pythia8: src/main_$$@.cpp $(PYTHIALIB) ../obj/libParticleEvent.so
	$(CXX) $^ -w $(INCLUDE) $(CXXFLAGS) $(ROOT) $(PYTHIA) -ldl -o $@.exe

test: src/main_$$@.cpp $(PYTHIALIB)
	$(CXX) $^ -w $(INCLUDE) $(CXXFLAGS) $(ROOT) $(PYTHIA) -ldl -o $@.exe

../obj/ParticleDct.cc: ../include/ParticleEvent.h ../include/ParticleLinkDef.h
	cd .. && $(MAKE) obj/ParticleDct.cc

# Objects

../obj/ParticleEvent.o: ../src/ParticleEvent.cpp
	cd .. && $(MAKE) obj/ParticleEvent.o

../obj/ParticleEventDict.o: ../obj/ParticleDct.cc
	cd .. && $(MAKE) obj/ParticleEventDict.o

# Shared objects

../obj/libParticleEvent.so: ../obj/ParticleEventDict.o ../obj/ParticleEvent.o
	cd .. && $(MAKE) obj/libParticleEvent.so

# Clean.
clean:
	rm -f *.log *.exe *.d *.so

