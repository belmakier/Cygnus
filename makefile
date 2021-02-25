CURDIR = $(shell pwd)
SRCDIR = $(CURDIR)/src
INCDIR = $(CURDIR)/include
BINDIR = $(CURDIR)/bin

ROOT_LIBS = `root-config --glibs` -lSpectrum -lTreePlayer -lMathMore

LIBRS = $(INCDIR) $(ROOT_LIBS) $(GSLLIBS)
INCLUDE = $(INCDIR)

CFLAGS = -std=c++11 -g -fPIC `root-config --cflags` `gsl-config --cflags` -Wno-unused-parameter

PLATFORM:=$(shell uname)
$(info PLATFORM: $(PLATFORM))
 
ifeq ($(PLATFORM),Darwin)
export __APPLE__:= 1
CFLAGS     += -Qunused-arguments
CPP        = clang++
else
export __LINUX__:= 1
CPP        = g++
endif

HEAD = $(wildcard include/*.h)
OBJECTS = $(patsubst include/%.h,lib/%.so,$(HEAD))

TARGET = bin/libCygnus.so

main: $(TARGET)
	@printf "Make complete\n"

$(TARGET): $(OBJECTS) lib/DictOutput.cxx lib bin
	@printf "Now compiling shared library $@\n"
	@$(CPP) $(CFLAGS) -I$(INCDIR) -I. -L$(LIBRS) -o $@ -shared lib/DictOutput.cxx $(OBJECTS) 

lib/DictOutput.cxx: $(HEAD)
	@printf "Linking libraries\n"
	@rootcint -f $@ -c -p $(HEAD) lib/linkdef.h

lib bin:
	@mkdir -p bin lib

lib/%.so: src/%.cxx include/%.h
	@printf "Now compiling library $@\n"
	@$(CPP) $(CFLAGS) -I$(INCDIR) -L$(LIBRS) -o $@ -shared -c $<
 
lib/%.so: src/*/%.cxx include/%.h	
	@printf "Now compiling library $@\n"
	@$(CPP) $(CFLAGS) -I$(INCDIR) -L$(LIBRS) -o $@ -share -c $< 

clean:  
	@printf "Tidying up...\n"
	@rm $(OBJECTS)
	@rm lib/DictOutput*
	@rm bin/*
	
