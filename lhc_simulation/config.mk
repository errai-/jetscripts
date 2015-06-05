# These may vary slightly depending on the system in use.

SHELL = /bin/sh

#######################################
# Using g++, the usual case:
#CXX := g++
#CFLAGS :=
# Using clang, more hardcore debugging:
CXX := clang++
CFLAGS :=
#######################################

CFLAGS := $(CFLAGS) -g -Wall -Wextra -Werror -Wshadow -Wno-error=unused-parameter -Wno-error=unused-function -Wno-error=unknown-pragmas -march=native -O3 -pedantic -m64
CXXFLAGS := $(CFLAGS) -std=c++11

PYTHIA6:=-lEG -lEGPythia6 -L/home/hannu/Cern/pythia6 -lPythia6
PYTHIA8:=$(shell pythia8-config --cflags --libs)
FASTJET:=$(shell fastjet-config --cxxflags --libs)
ROOT:=$(shell root-config --cflags --glibs) -lGenVector
ROOT_INCLUDE:=$(shell root-config --cflags)
CINTERPRET=$(ROOTSYS)/bin/rootcling
# ROOT 5:
#CINTERPRET=$(shell $(ROOTSYS)/bin/rootcint)

INCLUDE := -I./include
