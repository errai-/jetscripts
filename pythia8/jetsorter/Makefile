######################################################################
# Makefile for building Pythia-ROOT-Fastjet files
# Hannu Siikonen, 29.12.2014 
# Modified from pythia8 root examples makefile 
# (Axel Naumann, 2011-03-03)
######################################################################

# A few variables used in this Makefile:
SHELL = /bin/sh
CXXFLAGS = -O2 -ansi -pedantic -W -Wall -Wshadow -m64

EX           := hist tree
EXE          := $(addsuffix .exe,$(EX))
PYTHIA8      ?= /usr/local # valid if 'make install' has been performed 
STATICLIB    := $(PYTHIA8)/lib/archive/libpythia8.a
SHAREDLIB    := $(PYTHIA8)/lib/libpythia8.so
DICTCXXFLAGS := -I$(PYTHIA8)/include 
ROOTCXXFLAGS := $(DICTCXXFLAGS) $(shell root-config --cflags)

ROOTLIBS = $(shell root-config --libs)

# Libraries to include if GZIP support is enabled
ifeq (x$(ENABLEGZIP),xyes)
LIBGZIP=-L$(BOOSTLIBLOCATION) -lboost_iostreams -L$(ZLIBLOCATION) -lz
endif

# LDFLAGS1 for static library, LDFLAGS2 for shared library
LDFLAGS1 := $(shell root-config --ldflags --glibs) \
  -L$(PYTHIA8)/lib/archive -lpythia8 -llhapdfdummy $(LIBGZIP)
LDFLAGS2 := $(shell root-config --ldflags --glibs) \
  $(shell /home/hannu/pythia8/bin/pythia8-config --cxxflags --libs) \
  $(shell fastjet-config --cxxflags --libs) 

# Default target; make examples (but not shared dictionary)
all: $(EX)

jetsorter: $(SHAREDLIB)
	$(CXX) $(DICTCXXFLAGS)  $@.cc \
	$(CXXFLAGS) -Wno-shadow $(ROOTCXXFLAGS) \
		-o $@.exe $(LDFLAGS2)

# static PYTHIA 8 library

# Error messages if PYTHIA libraries don't exist
$(STATICLIB):
	@echo "Error: PYTHIA 8 archive library must be built"
	@false
$(SHAREDLIB):
	@echo "Error: PYTHIA 8 shared library must be built"
	@false

# Clean up
clean:
	rm -f $(EXE) *.log *.exe sortedjets.root *.d *.so \

