######################################################################
# Makefile for building Pythia's ROOT examples.
# Axel Naumann, 2011-03-03
######################################################################

# Need this to get SHAREDSUFFIX (e.g. dylib or so)
-include ../../config.mk

# A few variables used in this Makefile:
EX           := hist tree
EXE          := $(addsuffix .exe,$(EX))
PYTHIA8      ?= $(PWD)/..
STATICLIB    := $(PYTHIA8)/lib/archive/libpythia8.a
SHAREDLIB    := $(PYTHIA8)/lib/libpythia8.$(SHAREDSUFFIX)
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
  -L$(PYTHIA8)/lib -lpythia8 -llhapdfdummy $(LIBGZIP)
CXXFLAGS += -I/home/hannu/fastjet-install/include \
	-L$(PYTHIA8)/lib -lpythia8 -llhapdfdummy -L/home/hannu/fastjet-install/lib -lfastjet

# Default target; make examples (but not shared dictionary)
all: $(EX)

ptcut: $(STATICLIB)
	$(CXX) $(DICTCXXFLAGS)  $@.cc \
	`/home/hannu/fastjet-install/bin/fastjet-config --cxxflags --plugins` \
	$(CXXFLAGS) -Wno-shadow $(ROOTCXXFLAGS) \
		-o $@.exe $(LDFLAGS2)

# Rule to build tree example. Needs dictionary to be built and
# static PYTHIA 8 library
event2tree: $(STATICLIB)
	$(CXX) $(DICTCXXFLAGS)  $@.cc \
	`/home/hannu/fastjet-install/bin/fastjet-config --cxxflags --plugins` \
	$(CXXFLAGS) -Wno-shadow $(ROOTCXXFLAGS) \
		-o $@.exe $(LDFLAGS2)

ptcut_progressive: $(STATICLIB)
	$(CXX) $(DICTCXXFLAGS)  $@.cc \
	`/home/hannu/fastjet-install/bin/fastjet-config --cxxflags --plugins` \
	$(CXXFLAGS) -Wno-shadow $(ROOTCXXFLAGS) \
		-o $@.exe $(LDFLAGS2)

ptcut_stable: $(STATICLIB)
	$(CXX) $(DICTCXXFLAGS)  $@.cc \
	`/home/hannu/fastjet-install/bin/fastjet-config --cxxflags --plugins` \
	$(CXXFLAGS) -Wno-shadow $(ROOTCXXFLAGS) \
		-o $@.exe $(LDFLAGS2)

ptcut_stable_sum: $(STATICLIB)
	$(CXX) $(DICTCXXFLAGS)  $@.cc \
	`/home/hannu/fastjet-install/bin/fastjet-config --cxxflags --plugins` \
	$(CXXFLAGS) -Wno-shadow $(ROOTCXXFLAGS) \
		-o $@.exe $(LDFLAGS2)

ptcut_datahists: $(STATICLIB)
	$(CXX) $(DICTCXXFLAGS)  $@.cc \
	`/home/hannu/fastjet-install/bin/fastjet-config --cxxflags --plugins` \
	$(CXXFLAGS) -Wno-shadow $(ROOTCXXFLAGS) \
		-o $@.exe $(LDFLAGS2)

# Error messages if PYTHIA libraries don't exist
$(STATICLIB):
	@echo "Error: PYTHIA 8 archive library must be built"
	@false
$(SHAREDLIB):
	@echo "Error: PYTHIA 8 shared library must be built"
	@false

# Clean up
clean:
	rm -f $(EXE) hist.root pythiaDict.* \
               treeDict.cc treeDict.h pytree.root

