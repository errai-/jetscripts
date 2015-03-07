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

StoreParticles.so : src/StoreParticles.cc ../obj/libPrtclEvent.so
	$(CXX) -shared -fPIC $^ $(INCLUDE) $(ROOT) $(CXXFLAGS) $(LDFLAGS) -o $@  

../obj/PrtclDct.cc: ../include/PrtclEvent.h ../include/PrtclLinkDef.h
	cd .. && $(MAKE) obj/PrtclDct.cc

# Objects

../obj/PrtclEvent.o: ../src/PrtclEvent.cpp
	cd .. && $(MAKE) obj/PrtclEvent.o

../obj/PrtclEventDict.o: ../obj/PrtclDct.cc
	cd .. && $(MAKE) obj/PrtclEventDict.o

# Shared objects

../obj/libPrtclEvent.so: ../obj/PrtclEventDict.o ../obj/PrtclEvent.o
	cd .. && $(MAKE) obj/libPrtclEvent.so

clean:
	rm -f $(ALLCCFILES:.cc=.so)
