###################################################################
## Makefile for building ROOT-fastjet analysis of simulated data ## 
## Hannu Siikonen, 20.8.2015                                     ##  
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

all: pythia8 pythia6 libHerwigppTree.so jetanalysis ptcut

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

analysis: jetanalysis/main_standard.cpp jetbase.o analysis.o libPrtclEvent.so libJetEvent.so
	cd jetanalysis && $(CXX) main_standard.cpp jetbase.o analysis.o -L../lib -lPrtclEvent -lJetEvent -w \
	$(INCLUDE) $(CXXFLAGS) $(ROOT) $(FASTJET) -I. -ldl -o ../analysis.exe

jetanalysis.o: JetAnalysis.cpp libJetEvent.so libPrtclEvent.so
	cd jetanalysis && $(CXX) -c JetAnalysis.cpp -L../lib -lPrtclEvent -lJetEvent \
	$(INCLUDE) $(CXXFLAGS) -w $(ROOT) $(FASTJET) -ldl -o $@

analysis.o: Analysis.cpp jetbase.o
	cd jetanalysis && $(CXX) -c Analysis.cpp jetbase.o -L../lib -lPrtclEvent -lJetEvent \
	$(INCLUDE) $(CXXFLAGS) -w $(ROOT) $(FASTJET) -ldl -o $@

jetbase.o: JetBase.cpp libJetEvent.so libPrtclEvent.so
	cd jetanalysis && $(CXX) -c JetBase.cpp -L../lib -lPrtclEvent -lJetEvent \
	$(INCLUDE) $(CXXFLAGS) -w $(ROOT) $(FASTJET) -ldl -o $@

testjet: jetanalysis/main.cpp testanalysis.o libPrtclEvent.so libJetEvent.so
	cd jetanalysis && $(CXX) main.cpp testanalysis.o -L../lib -lPrtclEvent -lJetEvent -w \
	$(INCLUDE) $(CXXFLAGS) $(ROOT) $(FASTJET) -I. -ldl -o ../testjet.exe

testanalysis.o: JetTest.cpp libJetEvent.so libPrtclEvent.so
	cd jetanalysis && $(CXX) -c JetTest.cpp -L../lib -lPrtclEvent -lJetEvent \
	$(INCLUDE) $(CXXFLAGS) -w $(ROOT) $(FASTJET) -ldl -o $@


ptcut: ptcut/main.cpp ptcut.o libPrtclEvent.so
	cd ptcut && $(CXX) main.cpp ptcut.o -L../lib -lPrtclEvent -w \
	$(INCLUDE) $(CXXFLAGS) $(ROOT) $(FASTJET) -I. -ldl -o ../$@.exe

ptcut.o: PTCut.cpp libPrtclEvent.so
	cd ptcut && $(CXX) -c PTCut.cpp -L../lib -lPrtclEvent $(INCLUDE) \
	$(CXXFLAGS) -w $(ROOT) -ldl -o $@

# Pythia8

pythia8: pythia8/main.cpp pythia8.o libPrtclEvent.so
	cd pythia8 && $(CXX) main.cpp pythia8.o -w -L../lib -lPrtclEvent $(INCLUDE) \
	$(CXXFLAGS) $(ROOT) $(PYTHIA8) -ldl -o pythia8.exe

pythia8.o: Pythia8Tree.cpp libPrtclEvent.so
	cd pythia8 && $(CXX) -c Pythia8Tree.cpp -L../lib -lPrtclEvent \
	$(INCLUDE) $(CXXFLAGS) -w $(ROOT) $(PYTHIA8) -ldl -o $@

run_pythia8:
	bash bash/pythia8.sh $(NUM_EVT) $(MODE) $(NUM_PROC)

# Pythia6

pythia6: pythia6/main.cpp pythia6.o libPrtclEvent.so
	cd pythia6 && $(CXX) main.cpp pythia6.o -w -L../lib -lPrtclEvent $(INCLUDE) \
	$(CXXFLAGS) $(ROOT) $(LHAPDF) $(PYTHIA6) -ldl -o pythia6.exe

pythia6.o: Pythia6Tree.cpp libPrtclEvent.so
	cd pythia6 && $(CXX) -c Pythia6Tree.cpp -L../lib -lPrtclEvent \
	$(INCLUDE) $(CXXFLAGS) -w $(ROOT) $(PYTHIA6) -ldl -o $@

run_pythia6:
	bash bash/pythia6.sh $(NUM_EVT) $(MODE) $(NUM_PROC)

# Herwig

libHerwigppTree.so : HerwigppTree.cpp HerwigppTree.h libPrtclEvent.so
	cd herwig && $(CXX) -shared -fPIC HerwigppTree.cpp \
	-L../lib -lPrtclEvent -isystem $(shell herwig-config --includedir) $(INCLUDE) $(ROOT) $(CXXFLAGS) $(LDFLAGS) -o ../lib/$@

run_herwig:
	bash bash/herwig.sh $(NUM_EVT) $(MODE) $(NUM_PROC)

# Clean.
clean:
	rm -f *.log *.exe *.d *.so events/*.o jetsorting/*.o events/*Dct.* lib/* \
	herwig/*log herwig/*run herwig/*out herwig/*tex
