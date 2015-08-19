###################################################################
## Makefile for building ROOT-fastjet analysis of simulated data ## 
## Hannu Siikonen, 6.5.2015                                      ##  
###################################################################

include config.mk

.SECONDEXPANSION:
.DELETE_ON_ERROR:
.PHONY: all clean run_herwig pythia6 run_pythia6 pythia8 run_pythia8 jetanalysis ptcut
VPATH= events:jetanalysis:lib:pythia6:pythia8:herwig:ptcut

INCLUDE := $(INCLUDE) -I.

# Parameter vals:
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

jetanalysis: jetanalysis/main.cpp jetanalysis.o libPrtclEvent.so libJetEvent.so
	cd jetsorting && $(CXX) main.cpp jetanalysis.o -L../lib -lPrtclEvent -lJetEvent -w \
	$(INCLUDE) $(CXXFLAGS) $(ROOT) $(FASTJET) -I. -ldl -o ../analysis.exe

jetanalysis.o: JetAnalysis.cpp libJetEvent.so libPrtclEvent.so
	cd jetsorting && $(CXX) -c JetAnalysis.cpp -L../lib -lPrtclEvent -lJetEvent \
	$(INCLUDE) $(CXXFLAGS) -w $(ROOT) $(FASTJET) -ldl -o $@


ptcut: ptcut/main.cpp ptcut.o libPrtclEvent.so
	cd jetsorting && $(CXX) main.cpp ptcut.o -L../lib -lPrtclEvent -w \
	$(INCLUDE) $(CXXFLAGS) $(ROOT) $(FASTJET) -I. -ldl -o ../$@.exe

ptcut.o: PTCut.cpp libPrtclEvent.so
	cd jetsorting && $(CXX) -c PTCut.cpp -L../lib -lPrtclEvent $(INCLUDE) \
	$(CXXFLAGS) -w $(ROOT) -ldl -o $@

# Pythia8

pythia8: pythia8/main.cpp libPrtclEvent.so pythia8/functions.h
	cd pythia8 && $(CXX) main.cpp -w -L../lib -lPrtclEvent $(INCLUDE) \
	$(CXXFLAGS) $(ROOT) $(PYTHIA8) -ldl -o pythia8.exe

run_pythia8:
	bash bash/pythia8.sh $(NUM_EVT) $(MODE) $(NUM_PROC)

# Pythia6

pythia6: pythia6/main.cpp libPrtclEvent.so pythia6/functions.h
	cd pythia6 && $(CXX) main.cpp -w -L../lib -lPrtclEvent $(INCLUDE) \
	$(CXXFLAGS) $(ROOT) $(LHAPDF) $(PYTHIA6) -ldl -o pythia6.exe

run_pythia6:
	bash bash/pythia6.sh $(NUM_EVT) $(MODE) $(NUM_PROC)

# Herwig

libStoreParticles.so : StoreParticles.cc StoreParticles.h libPrtclEvent.so
	cd herwig && $(CXX) -shared -fPIC StoreParticles.cc \
	-L../lib -lPrtclEvent -isystem $(shell herwig-config --includedir) $(INCLUDE) $(ROOT) $(CXXFLAGS) $(LDFLAGS) -o ../lib/$@

run_herwig:
	bash bash/herwig.sh $(NUM_EVT) $(MODE) $(NUM_PROC)

# Clean.
clean:
	rm -f *.log *.exe *.d *.so events/*.o jetsorting/*.o events/*Dct.* lib/* \
	herwig/*log herwig/*run herwig/*out herwig/*tex
