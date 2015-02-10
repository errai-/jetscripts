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

StoreParticles.so : src/StoreParticles.cc ../obj/libP2TLib.so ../obj/minimalevent.so
	$(CXX) -shared -fPIC $^ $(INCLUDE) $(ROOT) $(CXXFLAGS) $(LDFLAGS) -o $@  

../obj/PTDct.cc:
	cd .. && $(MAKE) obj/PTDct.cc

../obj/libP2TLib.so: ../obj/PTDct.cc ../src/MinimalEvent.cpp 
	cd .. && $(MAKE) obj/libP2TLib.so

../obj/minimalevent.so: ../src/MinimalEvent.cpp 
	cd .. && $(MAKE) obj/minimalevent.so

clean:
	rm -f $(ALLCCFILES:.cc=.so)
