#############################################################
## Makefile for building Pythia-ROOT-Fastjet files (ptcut) ##
## Hannu Siikonen, 30.12.2014                              ##
#############################################################

# It is assumed that pythia8 and fastjet are installed under /usr/local
# or some other known location. If not, see Makefile_long for adding proper
# includes.

# Separate cxxflags are required only for root, 
# since root is not installed under /usr/local
CXXFLAGS := -O2 -ansi -pedantic -W -Wall -Wshadow -Wno-shadow -m64
ROOTCXXFLAGS := $(shell root-config --cflags)

# The libraries are loaded using proper config-files
SHELL = /bin/sh
LDFLAGS := $(shell root-config --ldflags --glibs) \
  $(shell pythia8-config --libs) $(shell fastjet-config --libs) 

ptcut: 
	$(CXX) $@.cc $(CXXFLAGS) $(ROOTCXXFLAGS) $(LDFLAGS) -o $@.exe 

ptcut_progressive: 
	$(CXX) $@.cc $(CXXFLAGS) $(ROOTCXXFLAGS) $(LDFLAGS) -o $@.exe 

ptcut_stable: 
	$(CXX) $@.cc $(CXXFLAGS) $(ROOTCXXFLAGS) $(LDFLAGS) -o $@.exe 

ptcut_stable_sum: 
	$(CXX) $@.cc $(CXXFLAGS) $(ROOTCXXFLAGS) $(LDFLAGS) -o $@.exe 

ptcut_datahists: 
	$(CXX) $@.cc $(CXXFLAGS) $(ROOTCXXFLAGS) $(LDFLAGS) -o $@.exe 

clean:
	rm -f $(EXE) *.log *.exe *.d *.so 

