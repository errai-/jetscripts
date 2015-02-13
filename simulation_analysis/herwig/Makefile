############################
# Herwig++ 2.7.1           #
# Hannu Siikonen 10.2.2015 #
############################

-include ../config.mk

ALLCCFILES=$(shell echo *.cc)

.SECONDEXPANSION:
.PHONY: all clean

all : $(ALLCCFILES:.cc=.so)

INCLUDE := $(INCLUDE) -I../include -I..

#%.so : %.cc %.h 
#	$(CXX) -shared -fPIC $(CPPFLAGS) $(ROOT) $(CXXFLAGS) $(LDFLAGS) $< -o $@  

# Main targets

StoreParticles.so : src/StoreParticles.cc ../obj/libSimEvent.so
	$(CXX) -shared -fPIC $^ $(INCLUDE) $(ROOT) $(CXXFLAGS) $(LDFLAGS) -o $@  

../obj/SimDct.cc: ../include/SimEvent.h ../include/SimLinkDef.h
	cd .. && $(MAKE) obj/SimDct.cc

# Objects

../obj/SimEvent.o: ../src/SimEvent.cpp
	cd .. && $(MAKE) obj/SimEvent.o

../obj/SimEventDict.o: ../obj/SimDct.cc
	cd .. && $(MAKE) obj/SimEventDict.o

# Shared objects

../obj/libSimEvent.so: ../obj/SimEventDict.o ../obj/SimEvent.o
	cd .. && $(MAKE) obj/libSimEvent.so

clean:
	rm -f $(ALLCCFILES:.cc=.so)
