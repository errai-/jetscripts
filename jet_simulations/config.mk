# These may vary slightly depending on the system in use.

SHELL = /bin/sh

#######################################
# Using g++, the usual case:
CXX = g++
CFLAGS :=
# Using clang, more hardcore debugging:
#CXX = clang
#CFLAGS := -v -fsanitize=undefined 
#######################################

CFLAGS := $(CFLAGS) -g -Wall -Wextra -Werror -Wshadow -Wno-error=unused-parameter -Wno-error=unused-function -Wno-error=unknown-pragmas -march=native -O3 -pedantic -m64
CXXFLAGS := $(CFLAGS) -std=c++11

ROOT=$(shell root-config --cflags --glibs) -lGenVector
PYTHIA6=-lEG -lEGPythia6 -L/home/hannu/Cern/pythia6 -lPythia6  
PYTHIA8=$(shell pythia8-config --cflags --libs)
FASTJET=$(shell fastjet-config --cxxflags --libs)

INCLUDE := -I./include