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



# 970821 Created by Max Moroz and Igor Markov,  VLSICAD ABKGroup, UCLA/CS
# 990812 ilm upgraded for use with g++

.SUFFIXES:        .o .cxx .c .exeS .exe .h.gch

# ##################
# standard variables
# ##################

 OPTFLAGS        = -g -DABKDEBUG 
#OPTFLAGS        = -O3 # use this for optimization instead of the above
CCFLAGS          = -fPIC
DEPENDOPT        = -xM1
CC               = g++ $(OPTFLAGS) $(CCFLAGS)
C                = g++
LD               = g++ -shared
LDFINAL          = g++ -Wl,--no-as-needed
CC4LIBS          = CC6
AR               = ar cr
SODIRKEY =  -Wl,-rpath,
SONAMEKEY        = -Wl,-soname,
# change to -h if you need "immobile" dynamic libraries.
                      # LD_LIBRARY_PATH overrides -h
                      # Make sure there is no conflict

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

# #########################################################
# standard variables that need modification in each project
# #########################################################

THISLIB = MLPart
INCLUDES    = -I$(UCLA)  $(STL)

OBJS             = mlParams.o baseMLPart.o doML.o mlPart.o callPartitioner.o

SRCS        =  $(OBJS:%.o=%.cxx)
LIBPATHS    = $(SODIRKEY)$(LIBBASE)

LIBS             =-L$(LIBBASE) $(SODIRKEY). $(LIBPATHS) -l$(THISLIB)\
		   -Wl,--start-group -lABKCommon -lCtainers\
                   -lCombi -lStats -lPlacement -lHGraph -lGeoms\
                   -lPartEvals -lPartLegality -lPartitioning\
                   -lPartitioners -lClusteredHGraph -lFMPart\
                   -lSmallPart -lGenHier -lFilledHier -Wl,--end-group

SLIBS            =-L$(SLIBBASE) -l$(THISLIB) \
                  -Wl,--start-group -lFMPart -lSmallPart -lPartitioners -lPartEvals\
                  -lPartLegality -lClusteredHGraph -lFilledHier -lGenHier\
                  -lPartitioning -lHGraph -lGeoms -lPlacement -lStats\
                  -lCtainers -lCombi -lABKCommon -Wl,--end-group

PURIFY  = purify -cache_dir=cache -always_use_cache_dir\
	  -windows=no -log_file=purify.log

# ################
# standard targets
# ################

all:    lib

lib: $(ARNAME) $(SOFULL)

$(SOFULL):        $(OBJS)
	/bin/rm -f $(SOFULL)
	$(LD) $(INCLUDES) -o $(SOFULL) $(OBJS) $(SOLABEL)

$(ARNAME):        $(OBJS)
	/bin/rm -f $(ARNAME)
	$(AR) $(ARNAME) $(OBJS) $(RNGOBJS)    

.cxx.o:
	$(CC) $(INCLUDES) -c $*.cxx

.c.o:
	$(C) $(INCLUDES) -c $*.c

.h.h.gch:
	$(CC) $(INCLUDES) $(ABKSYMBOLS) $*.h

%.exe: %.cxx %.o
	/bin/rm -f $*.exe
	$(LDFINAL) -o $*.exe $*.o $(LIBS)

%.exeS: %.cxx %.o
	/bin/rm -f $*.exeS
	$(LDFINAL) -o $*.exeS $*.o -static $(SLIBS)


depend:
	/bin/rm -f depend.mk

depend.mk:
	$(CC) $(DEPENDOPT) $(INCLUDES) $(SRCS) *.cxx >depend.mk

clean:        
	rm -rf *.o Templates.DB *.exe* $(SOFULL) $(ARNAME) core ir.out\
	lib*pure* cache SunWS_cache *notime new.out *.h.gch

TEST_EXES = MLPartTest0.exe  MLPartTest1.exe  MLPartTest2.exe  MLPartTest6.exe\
       MLPartTest7.exe MLPartTest8.exe
test: lib $(TEST_EXES)

# Test executables link against this package's library; without this,
# "make -jN test" can link a test before "lib" has (re)built it.
$(TEST_EXES): lib

testS: lib MLPartTest0.exeS  MLPartTest1.exeS  MLPartTest2.exeS  MLPartTest6.exeS\
       MLPartTest7.exeS MLPartTest8.exeS

# include depend.mk
