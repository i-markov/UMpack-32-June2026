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

# 971006 ilm added dependency on Stats
# 990812 ilm upgraded for use with g++  

.SUFFIXES:        .o .cxx .c .exe .exeS .h.gch

# ##################
# standard variables
# ##################

 OPTFLAGS        = -g -DABKDEBUG -DHGDEBUG -DPARTDEBUG -DCAPODEBUG $(ABKSYMBOLS)
#OPTFLAGS        = -O3 # use this for optimization instead of the above
CCFLAGS          = -fPIC
DEPENDOPT        = -xM1
CC               = g++ $(OPTFLAGS) $(CCFLAGS)
C                = gcc $(OPTFLAGS) $(CCFLAGS)
LD               = g++ -shared
LDFINAL          = g++ -Wl,--no-as-needed
CC4LIBS          = g++
# could be g++
AR               = ar cr
SODIRKEY         = -Wl,-rpath,
SONAMEKEY        = -Wl,-soname,
# change to -h if you need "immobile"dynamic libraries.
# LD_LIBRARY_PATH overrides -h
# Make sure there is no conflict
             
UCLA     := $(abspath $(dir $(abspath $(lastword $(MAKEFILE_LIST))))/..)
IMPORT   = NOT_USED
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

NSLLIBS     =  -lnsl -lsocket -lintl


# #########################################################
# standard variables that need modification in each project
# #########################################################

THISLIB = Capo
INCLUDES= $(STL) -I$(UCLA)

QLIBS  = 
QLIBS  =   
YACC_LIB         = 
#

LIBPATHS         = $(SODIRKEY)$(LIBBASE)

LIBS             =-L$(LIBBASE) $(SODIRKEY). $(LIBPATHS) \
                  -Wl,--start-group -l$(THISLIB) \
                  -lMLPart -lFMPart -lPartEvals -lPartitioning -lPartitioners \
                  -lPartLegality $(QLIBS) \
                  $(YACC_LIB) -lConstraints -lPlacement -lStats -lHGraph \
                  -lCombi -lPlaceEvals -lClusteredHGraph -lSmallPart \
                  -lSmallPlacers -lRBPlace -lHGraphWDims \
                  -lCongestionMaps -lParquetDBFromRBP -lParquetFP \
                  -lGeomTrees -lGeoms -lCtainers -lSmallPlacement \
                  -lShellPart -lGenHier -lFilledHier -lAnalytPl \
                  -lABKCommon -Wl,--end-group

# -lSTA -lTDplace -lBoxPlace -lLineSearch -lLinAlg \

SLIBS            =-L$(SLIBBASE) \
                  -Wl,--start-group -l$(THISLIB) -lShellPart \
		  -lSmallPlacers -lSmallPlacement -lRBPlace -lMLPart \
		  -lFMPart -lSmallPart -lPartitioners \
		  -lClusteredHGraph -lGenHier -lFilledHier -lPartEvals \
		  -lPartLegality $(QLIBS) \
		  -lPlaceEvals -lPartitioning -lHGraphWDims -lHGraph \
		  $(YACC_LIB) -lParquetDBFromRBP -lParquetFP \
		  -lConstraints -lPlacement -lGeomTrees -lGeoms -lStats \
		  -lCtainers -lCombi -lStats -lAnalytPl -lCongestionMaps \
                  -lABKCommon -Wl,--end-group

# -lSTA -lTDplace -lBoxPlace -lLineSearch -lLinAlg \
# -lefence
 
OBJS	 = capoParams.o capoBinPl.o baseBinSplitter.o \
	   splitSmallH.o splitSmallV.o\
	   splitLargeH.o splitLargeV.o\
	   partProbForCapo.o partProbCtor.o \
	   capoBin.o capoPlacer.o capoPrinting.o\
  	   constrForCapo.o refineBin.o \
	   capoMainLoop.o overlapping.o weber.o\
	   subHGForCapo.o laSplitter.o splitRow.o\
	   constrainingHG.o thunderPart.o cutPosition.o\
	   cutShifting.o

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
 
test:  CapoTest0.exe CapoTest1.exe
testS: CapoTest1.exeS
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
#TEST_EXP       = aztec.out.D.A.dag
TEST_DIFF       = test.diff
#TEST_ARGS      = -aux aztec.aux -wm D -fm A
TEST_ARGS       = -f mini2


# ########################
# project-specific targets
# ########################

# Using a C-shell regression script is recommended instead
# (see file "regression")

regression:        CapoTest0.exe CapoTest1.exe
	./regression

correct:
	$(TESTNAME) $(TEST_ARGS) > $(TEST_EXP) < $(TEST_IN)
	perl -n -i -e "print unless /Created|sec/;" $(TEST_EXP)
