# These may vary slightly depending on the system in use.

SHELL = /bin/sh

#######################################
# Using g++, the usual case:
CXX := g++
CFLAGS :=
# Using clang, more hardcore debugging:
#CXX := clang++
#CFLAGS :=
#######################################
# -Wno-error=unused-parameter -Wno-error=unused-function -Wno-error=unknown-pragmas
CFLAGS := $(CFLAGS) -g -Wall -Wextra -Werror -Wshadow -pedantic -m64
CXXFLAGS := $(CFLAGS) -std=c++11 -march=native -O2

PYTHIA6:=-lEG -lEGPythia6 -L/afs/cern.ch/user/h/hsiikone/Cern/pythia6 -lPythia6
//PYTHIA6:=-lEG -lEGPythia6 -L/home/hannu/Cern/pythia6 -lPythia6
PYTHIA8:=$(shell pythia8-config --cflags --libs)
FASTJET:=$(shell fastjet-config --cxxflags --libs)
ROOT:=$(shell root-config --cflags --glibs) -lGenVector
ROOT_INCLUDE:=$(shell root-config --cflags)
ROOT_SAFE:=$(shell root-config --glibs) -lGenVector -I$(shell root-config --incdir)
CINTERPRET=$(ROOTSYS)/bin/rootcling
# ROOT 5:
#CINTERPRET=$(shell $(ROOTSYS)/bin/rootcint)
LHAPDF:=$(shell lhapdf-config --cflags --libs)

INCLUDE := -I./include
