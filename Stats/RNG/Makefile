#/**************************************************************************
#***    
#*** Copyright (c) 1995-2000 Regents of the University of California,
#***               Andrew E. Caldwell, Andrew B. Kahng and Igor L. Markov
#*** Copyright (c) 2000-2012 Regents of the University of Michigan,
#***               Saurabh N. Adya, Jarrod A. Roy, David A. Papa and
#***               Igor L. Markov
#***
#***  Contact author(s): abk@cs.ucsd.edu, imarkov@umich.edu
#***  Original Affiliation:   UCLA, Computer Science Department,
#***                          Los Angeles, CA 90095-1596 USA
#***
#***  Permission is hereby granted, free of charge, to any person obtaining 
#***  a copy of this software and associated documentation files (the
#***  "Software"), to deal in the Software without restriction, including
#***  without limitation 
#***  the rights to use, copy, modify, merge, publish, distribute, sublicense, 
#***  and/or sell copies of the Software, and to permit persons to whom the 
#***  Software is furnished to do so, subject to the following conditions:
#***
#***  The above copyright notice and this permission notice shall be included
#***  in all copies or substantial portions of the Software.
#***
#*** THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, 
#*** EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
#*** OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. 
#*** IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
#*** CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT
#*** OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
#*** THE USE OR OTHER DEALINGS IN THE SOFTWARE.
#***
#***
#***************************************************************************/


# 000430 Created by Sachin Bansal,  VLSICAD ABKGroup, UCLA/CS

.SUFFIXES:       .a .o .cxx .c .cpp


# ##################
# standard variables
# ##################

 OPTFLAGS        = -g -DABKDEBUG 
#OPTFLAGS        = -O3 # use this for optimization instead of the above
CCFLAGS          = -fPIC
DEPENDOPT        = -xM1
CC               = g++ $(OPTFLAGS) $(CCFLAGS)
LD               = g++ -shared
LDFINAL          = g++ -Wl,--no-as-needed
CC4LIBS          = CC5
AR               = ar cr
SODIRKEY =  -Wl,-rpath,
SONAMEKEY        = -Wl,-soname,

UCLA     := $(abspath $(dir $(abspath $(lastword $(MAKEFILE_LIST))))/../..)
IMPORT   = NOT_USED
LIBBASE  := $(UCLA)/lib
SLIBBASE    = $(LIBBASE)
ARNAME      = lib$(THISLIB).a
SONAME      = lib$(THISLIB).so
SODIR       = $(CURDIR)
SOFULL      = $(SODIR)/$(SONAME)
SOLABEL     = $(SONAMEKEY)$(SONAME)
TESTNAME    = $(THISLIB)Test

STL         = 
# -I$(IMPORT)/STL/CC
COMMON      = $(UCLA)/ABKCommon

# #########################################################
# standard variables that need modification in each project
# #########################################################

THISLIB = Rand
INCLUDES= -I$(COMMON) 
OBJS    =  hist.o newran.o extreal.o myexcept.o
SRCS                =  $(OBJS:%.o=%.cpp)
LIBPATHS    = $(SODIRKEY)$(LIBBASE)
LIBS        = -L. -l$(THISLIB) $(LIBPATHS) -L$(LIBBASE)  -lcommon 
SLIBS       = $(ARNAME)                    -L$(SLIBBASE) -lcommon 
TEST_EXE    = rngsmoke.exe
TEST_OBJ    = rngsmoke.o

PURIFY  = purify -cache_dir=cache -always_use_cache_dir\
	  -windows=no -log_file=purify.log

# ################
# pakcage specific variables
# ################

OBJ_TR = tryrand.o tryrand1.o tryrand2.o tryrand3.o tryrand4.o hist.o \
      newran.o extreal.o myexcept.o

newranxx = include.h newran.h boolean.h myexcept.h


# ################
# standard targets
# ################

all:        lib

lib: $(SOFULL) $(ARNAME)

$(SOFULL):        $(OBJS)
	/bin/rm -f $(SOFULL)
	$(LD) $(INCLUDES) -o $(SOFULL) $(OBJS) $(SOLABEL)

$(ARNAME):        $(OBJS)
	/bin/rm -f $(ARNAME)
	$(AR) $(ARNAME) $(OBJS)     

%.o: %.cxx
	$(CC) $(INCLUDES) -c $<

depend:
	/bin/rm -f depend.mk

depend.mk:
	$(CC) $(DEPENDOPT) $(INCLUDES) $(SRCS) >depend.mk

clean:        
	rm -rf *.o Templates.DB $(TESTNAME)* $(TEST_EXE) $(SOFULL) $(SONAME) $(ARNAME)\
	core ir.out lib*pure* plot?.xy cache *.a tryrand SunWS_cache

test: lib $(TEST_EXE)
	./$(TEST_EXE) >/dev/null

# Test executables link against this package's library; without this,
# "make -jN test" can link a test before "lib" has (re)built it.
$(TEST_EXE): lib

testS: lib $(TEST_EXE)

tryrand: $(OBJ_TR) $(newranxx)
	$(LDFINAL) -o $@ $(OBJ_TR) -L. -lm

$(TEST_EXE): $(TEST_OBJ) $(OBJS)
	$(LDFINAL) -o $@ $< $(OBJS) -L. -lm

$(TEST_OBJ): rngsmoke.cxx
	$(CC) $(INCLUDES) -c $<

# include depend.mk
