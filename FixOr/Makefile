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

# 971006 ilm added dependency on StazenPlotz
# 971110 mro adapted for separation of ParserLEFDEF project
# 990812 ilm upgraded for use with g++

.SUFFIXES:        .o .cxx .c .exe .exeS .h.gch

# ##################
# standard variables
# ##################

 OPTFLAGS        = -g -DABKDEBUG $(ABKSYMBOLS)
#OPTFLAGS        = -O3 # use this for optimization instead of the above
CCFLAGS          = -fPIC
DEPENDOPT        = -MM
CC               = g++ $(OPTFLAGS) $(CCFLAGS)
C                = gcc $(OPTFLAGS) $(CCFLAGS)
LD               = g++ -shared
LDFINAL          = g++ -Wl,--no-as-needed
CC4LIBS          = g++
AR               = ar cr
SODIRKEY         = -Wl,-rpath,
SONAMEKEY        = -Wl,-soname,
# change to -h if you need "immobile" dynamic libraries.
# LD_LIBRARY_PATH overrides -h
# Make sure there is no conflict

UCLA     := $(abspath $(dir $(abspath $(lastword $(MAKEFILE_LIST))))/..)
LIBBASE  := $(UCLA)/lib
SLIBBASE         = $(LIBBASE)
STL              =
# -I$(IMPORT)/STL/CC

ARNAME           = lib$(THISLIB).a
SONAME           = lib$(THISLIB).so
SODIR            = $(CURDIR)
SOFULL           = $(SODIR)/$(SONAME)
SOLABEL          = $(SONAMEKEY)$(SONAME)
TESTNAME         = $(THISLIB)Test

# #########################################################
# standard variables that need modification in each project
# #########################################################

THISLIB = FixOr
INCLUDES= $(STL) -I$(UCLA)

QLIBS  =

YACC_LIB         =
#

LIBPATHS         = $(SODIRKEY)$(LIBBASE)

LIBS             =-L$(LIBBASE) $(SODIRKEY). $(LIBPATHS) \
                  -Wl,--start-group -l$(THISLIB) \
                  -lHGraphWDims -lConstraints -lPlacement -lStats -lHGraph \
                  $(QLIBS) $(YACC_LIB) -lCombi -lGeomTrees -lGeoms -lCtainers \
                  -lRBPlace -lCongestionMaps -lPlaceEvals \
                  -lABKCommon -Wl,--end-group

SLIBS            =-L$(SLIBBASE) \
                  -Wl,--start-group -l$(THISLIB) \
                  -lHGraphWDims -lConstraints -lPlacement -lStats -lHGraph \
                  $(QLIBS) $(YACC_LIB) -lCombi -lGeomTrees -lGeoms -lCtainers \
                  -lRBPlace -lCongestionMaps -lPlaceEvals \
                  -lABKCommon -Wl,--end-group

OBJS	 = greedyOrientOpt.o orientOptXFace.o

SRCS             =  $(OBJS:%.o=%.cxx)

PURIFY  = purify -cache_dir=cache -always_use_cache_dir\
	  -windows=no -log_file=purify.log

# The four definitions below are for use with the PRINT_VERSION_INFO macro
# defined in ABKCommon/abkversion.h
ABKPACKAGE =`pwd | sed 's/.*\///g'`
ABKSYSTEM  ="`uname -n`\040running\040`uname -sr | sed 's/ /-/g'`"
ABKUSER    =`whoami`@`uname -n`
ABKSYMBOLS =-D__ABKPACKAGE__=\"$(ABKPACKAGE)\" -D__ABKSYSTEM__=\"$(ABKSYSTEM)\"\
        -D__ABKUSER__=\"$(ABKUSER)\" -D__ABKLIBBASE__=\"$(LIBBASE)\"
# add -(ABKSYMBOLS) to the CC call defined by the .cxx.o rule

# ################
# standard targets
# ################

all:        lib

lib: $(ARNAME) $(SOFULL)

$(SOFULL):        $(OBJS)
	/bin/rm -f $(SOFULL)
	$(LD) $(INCLUDES) -o $(SOFULL) $(OBJS) $(SOLABEL)

$(ARNAME):        $(OBJS)
	/bin/rm -f $(ARNAME)
	$(AR) $(ARNAME) $(OBJS)

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
	$(CC) $(DEPENDOPT) $(INCLUDES) $(SRCS) >depend.mk

clean:
	rm -rf *.o Templates.DB *.exe* $(SOFULL) $(ARNAME) core ir.out\
	$(TEST_OUT) $(TEST_DIFF) lib*pure* cache SunWS_cache new.out *notime *~\
	*.h.gch

test:  FixOrTest1.exe
testS: FixOrTest1.exeS
.PHONY: depend, new, clean, test, lib, regression

# Uncomment for automatic handling of dependencies
# (use make depend to [re-]build dependencies)

# include depend.mk

# #################################
# sample project-specific variables
# #################################

TEST_OUT        = test.out
TEST_IN         = test.in
TEST_EXP        = test.exp
TEST_DIFF       = test.diff
TEST_ARGS       = -f mini2

# ########################
# project-specific targets
# ########################

# Using a shell regression script is recommended instead
# (see file "regression")

regression:        FixOrTest1.exe
	./regression

correct:
	$(TESTNAME) $(TEST_ARGS) > $(TEST_EXP) < $(TEST_IN)
	perl -n -i -e "print unless /Created|sec/;" $(TEST_EXP)
