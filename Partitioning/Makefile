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



.SUFFIXES:  .o .cxx  .exe .exeS .h.gch

# ##################
# standard variables
# ##################

# 990812 ilm upgraded for use with g++

OPTFLAGS    = -DABKDEBUG -g
#OPTFLAGS    = -O3
CCFLAGS     = -fPIC
DEPENDOPT   = -xM1
CC          = g++ $(OPTFLAGS) $(CCFLAGS)
LD          = g++ -shared
LDFINAL     = g++ -Wl,--no-as-needed
CC4LIBS          = CC6
AR          = ar cr
SODIRKEY =  -Wl,-rpath,
SONAMEKEY        = -Wl,-soname,
        
UCLA     := $(abspath $(dir $(abspath $(lastword $(MAKEFILE_LIST))))/..)
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

# #########################################################
# standard variables that need modification in each project
# #########################################################

THISLIB = Partitioning
OBJS        =  partProb.o partitionData.o partProbIO.o termiProp.o \
              relaxPart.o  subPartProb.o
SRCS        =  $(OBJS:%.o=%.cxx)
INCLUDES    = -I$(UCLA) -I. $(STL)
LIBPATHS    = $(SODIRKEY)$(LIBBASE)
LIBS        = -L. $(SODIRKEY). -l$(THISLIB) $(LIBPATHS) -L$(LIBBASE)\
	      -lABKCommon -lCtainers -lHGraph -lStats -lPlacement -lCombi -lGeoms
SLIBS       = $(ARNAME)                    -L$(SLIBBASE)\
	      -lABKCommon -lCtainers -lHGraph -lStats -lPlacement -lCombi -lGeoms
            
PURIFY      = purify -cache_dir=cache -always_use_cache_dir -windows=no\
              -log_file=purify.log 
# -append_logfile

# ################
# standard targets
# ################

all: lib

lib: $(ARNAME) $(SOFULL)

$(SOFULL):        $(OBJS)
	/bin/rm -f $(SOFULL)
	$(LD) $(INCLUDES) -o $(SOFULL) $(OBJS) $(SOLABEL)

$(ARNAME):        $(OBJS)
	/bin/rm -f $(ARNAME)
	$(AR) $(ARNAME) $(OBJS)

.cxx.o: 
	$(CC) $(INCLUDES) -c $*.cxx

.h.h.gch:
	$(CC) $(INCLUDES) $(ABKSYMBOLS) $*.h

%.exe: %.cxx %.o
	/bin/rm -f $@
	$(LDFINAL) -o $@ $*.o $(LIBS)

%.exeS: %.cxx %.o
	/bin/rm -f $@
	$(LDFINAL) -o $@ $*.o $(SLIBS)


depend:
	/bin/rm -f depend.mk

depend.mk:
	$(CC) $(DEPENDOPT) $(INCLUDES) $(SRCS) *.cxx >depend.mk

clean:  
	rm -rf *.o Templates.DB seeds.{in,out} *.exe* *.diff\
	$(SOFULL) $(SONAME) $(ARNAME) core ir.out *.notime depend.mk\
        cache SunWS_cache *notime new.out *.h.gch

TEST_EXES = PartitioningTest0.exe PartitioningTest1.exe PartitioningTest2.exe \
      PartitioningTest3.exe PartitioningTest4.exe PartitioningTest5.exe \
      PartitioningTest6.exe PartitioningTest7.exe PartitioningTest8.exe \
      PartitioningTest9.exe PartitioningTest10.exe PartitioningTest11.exe \
      PartitioningTest12.exe PartitioningTest13.exe
test: lib $(TEST_EXES)

# Test executables link against this package's library; without this,
# "make -jN test" can link a test before "lib" has (re)built it.
$(TEST_EXES): lib

testS: lib PartitioningTest0.exeS PartitioningTest1.exeS

.PHONY: depend, new, clean, test, lib

#include depend.mk

# ##################################################
# special variables related only the current project
# ##################################################

# ########################
# project-specific targets
# ########################
