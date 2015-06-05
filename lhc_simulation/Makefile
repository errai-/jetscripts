###################################################################
## Makefile for building ROOT-fastjet analysis of simulated data ## 
## Hannu Siikonen, 6.5.2015                                      ##  
###################################################################

include config.mk

.SECONDEXPANSION:
.DELETE_ON_ERROR:
.PHONY: all clean herwig run_herwig pythia6 run_pythia6 pythia8 run_pythia8 jetanalysis ptcut
VPATH= events:jetsorting:lib:pythia6:pythia8:herwig

INCLUDE := $(INCLUDE) -I.

# Parameter vals:

# The LHC-UE7-2 setting file, used by ATLAS
HFILE := LHC-UE7

NUM_EVT := 1000
MODE := 0
NUM_PROC := 1

all:
	@echo "No make for all objects provided."

# Generic

## Particle

libPrtclEvent.so: PrtclEvent.o PrtclEventDict.o
	cd events && $(CXX) -fPIC -shared $(CXXFLAGS) PrtclEvent.o \
	PrtclEventDict.o -o ../lib/$@

PrtclEvent.o: PrtclEvent.cpp PrtclEvent.h
	cd events && $(CXX) PrtclEvent.cpp -fPIC -pthread $(INCLUDE) \
	$(CXXFLAGS) $(ROOT) -c -w -ldl -o $@

PrtclEventDict.o: PrtclDct.cc
	cd events && $(CXX) $(INCLUDE) $(CXXFLAGS) $(ROOT_INCLUDE) -fPIC -pthread \
	-c ../events/PrtclDct.cc -o $@

PrtclDct.cc: PrtclEvent.h PrtclLinkDef.h
	cd events && $(CINTERPRET) -f $@ -c -p $(INCLUDE) $(CFLAGS) \
	PrtclEvent.h PrtclLinkDef.h && mv PrtclDct_rdict.pcm ../lib/.

## Jet

libJetEvent.so: JetEvent.o JetEventDict.o
	cd events && $(CXX) -fPIC -shared $(CXXFLAGS) JetEvent.o JetEventDict.o \
	-o ../lib/$@

JetEvent.o: JetEvent.cpp JetEvent.h
	cd events && $(CXX) -fPIC -pthread JetEvent.cpp $(INCLUDE) $(CXXFLAGS) \
	$(ROOT) -c -w -ldl -o $@

JetEventDict.o: JetDct.cc
	cd events && $(CXX) $(INCLUDE) $(ROOT_INCLUDE) $(CXXFLAGS) -fPIC -pthread \
	-c ../events/JetDct.cc -o $@

JetDct.cc: JetEvent.h JetLinkDef.h
	cd events && $(ROOTSYS)/bin/rootcint -f $@ -c -p $(INCLUDE) $(CFLAGS) \
	JetEvent.h JetLinkDef.h && mv JetDct_rdict.pcm ../lib/.

# Jet analysis

jetanalysis: jetsorting/main_$$@.cpp jetanalysis.o libPrtclEvent.so libJetEvent.so
	cd jetsorting && $(CXX) main_$@.cpp jetanalysis.o -L../lib -lPrtclEvent -lJetEvent -w \
	$(INCLUDE) $(CXXFLAGS) $(ROOT) $(FASTJET) -I. -ldl -o ../jetanalysis.exe

jetanalysis.o: jetsorting/JetAnalysis.cpp libJetEvent.so libPrtclEvent.so
	cd jetsorting && $(CXX) -c JetAnalysis.cpp -L../lib -lPrtclEvent -lJetEvent \
	$(INCLUDE) $(CXXFLAGS) -w $(ROOT) -ldl -o $@


ptcut: jetsorting/main_$$@.cpp ptcut.o libPrtclEvent.so
	cd jetsorting && $(CXX) main_$@.cpp ptcut.o -L../lib -lPrtclEvent -w \
	$(INCLUDE) $(CXXFLAGS) $(ROOT) $(FASTJET) -I. -ldl -o ../$@.exe

ptcut.o: jetsorting/PTCut.cpp libPrtclEvent.so
	cd jetsorting && $(CXX) -c PTCut.cpp -L../lib -lPrtclEvent $(INCLUDE) \
	$(CXXFLAGS) -w $(ROOT) -ldl -o $@

# Pythia8

pythia8: main_$$@.cpp libPrtclEvent.so pythia8_functions.h
	cd pythia8 && $(CXX) main_$@.cpp -w -L../lib -lPrtclEvent $(INCLUDE) \
	$(CXXFLAGS) $(ROOT) $(PYTHIA8) -ldl -o pythia8.exe

run_pythia8:
	bash pythia8.sh $(NUM_EVT) $(MODE) $(NUM_PROC)

# Pythia6

pythia6: main_$$@.cpp libPrtclEvent.so pythia6_functions.h
	cd pythia6 && $(CXX) main_$@.cpp -w -L../lib -lPrtclEvent $(INCLUDE) \
	$(CXXFLAGS) $(ROOT) $(PYTHIA6) -ldl -o pythia6.exe

run_pythia6:
	bash pythia6.sh $(NUM_EVT) $(MODE) $(NUM_PROC)

# Herwig

herwig: libStoreParticles.so
	cd herwig && Herwig++ read $(HFILE).in

run_herwig:
	bash herwig.sh $(NUM_EVT) $(MODE) $(NUM_PROC) $(HFILE)

libStoreParticles.so : StoreParticles.cc libPrtclEvent.so
	cd herwig && $(CXX) -shared -fPIC StoreParticles.cc \
	-L../lib -lPrtclEvent $(INCLUDE) $(ROOT) $(CXXFLAGS) $(LDFLAGS) -o ../lib/$@


# Clean.
clean:
	rm -f *.log *.exe *.d *.so events/*.o jetsorting/*.o events/*Dct.* lib/* \
	herwig/*log herwig/*run herwig/*out herwig/*tex
