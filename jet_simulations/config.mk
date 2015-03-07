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

CFLAGS := $(CFLAGS) -O2 -ansi -pedantic -W -Wall -Wshadow -m64 
CXXFLAGS := $(CFLAGS) -std=c++11

ROOT=$(shell root-config --cflags --glibs)
PYTHIA=$(shell pythia8-config --cflags --libs)
FASTJET=$(shell fastjet-config --cxxflags --libs)

PYTHIALIB=/usr/local/lib/libpythia8.so

INCLUDE := -I./include