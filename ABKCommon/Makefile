#/**************************************************************************
#***    
#*** Copyright (c) 1995-2000 Regents of the University of California,
#***               Andrew E. Caldwell, Andrew B. Kahng and Igor L. Markov
#*** Copyright (c) 2000-2012 Regents of the University of Michigan,
#***               Saurabh N. Adya, Jarrod A. Roy, David A. Papa and
#***               Igor L. Markov
#***
#***  Contact author(s): abk@cs.ucsd.edu, igor.markov1@gmail.com
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


# 970801 Created and tained by Igor Markov,  VLSICAD ABKGroup, UCLA/CS
# 970814 ilm added ir.out to clean, added +p switch, changed DEBUG to ABKDEBUG
# 970827 ilm added abkmanip.o to OBJS
# 970827 ilm added abkrand.o to OBJS
# 971113 mro added seed.o and old.o to OBJS (abkrand.o has been removed)
# 971114 mro added test9 to test use of seeds.in
# 971114 mro removed old.o from OBJS
# 971116 ilm added verbosity.o to OBJS
# 971116 ilm added test10
# 971216 mro added static linking
# 980212 ilm renamed abkmanip.o into abkio.o
# 980708 mro added md5our.o and test13
# 990208 ilm added definitions of symbols in $ABKSYMBOLS 
#            for use with PRINT_VERSION_INFO in abkversion.h
# 990812 ilm upgraded for use with g++

.SUFFIXES:	.o .cxx .exe .exeS .h.gch

 OPTFLAGS = -g -DABKDEBUG
#OPTFLAGS= -O3 # use this for optimization instead of the above
CCFLAGS  = -fPIC
# +w2 
# comment +w/+w2 when compiling with optimization
DEPENDOPT= -xM1
CC       = g++ $(OPTFLAGS) $(CCFLAGS)
LD       = g++ -shared
LDFINAL  = g++ -Wl,--no-as-needed
AR       = ar cr
SODIRKEY =  -Wl,-rpath,
SONAMEKEY= -Wl,-soname,
# change to -h if you need "immobile" dynamic libraries.
              # LD_LIBRARY_PATH overrides -h
              # Make sure there is no conflict

PURIFY   = purify -cache_dir=cache -windows=no -log_file=purify.log\
	   -always_use_cache_dir

UCLA     := $(abspath $(dir $(abspath $(lastword $(MAKEFILE_LIST))))/..)
IMPORT   = NOT_USED
CC4LIBS  = g++
LIBBASE  := $(UCLA)/lib
SLIBBASE = $(LIBBASE)
STL      = 
# -I$(IMPORT)/STL/CC
THISLIB = ABKCommon
ARNAME   = lib$(THISLIB).a
SONAME   = lib$(THISLIB).so
SODIR    = $(CURDIR)
TESTNAME = $(THISLIB)Test
MD5DIR   = $(CURDIR)/MD5


INCLUDES = $(STL) -I. -I$(UCLA)
LIBPATHS = $(SODIRKEY)$(CURDIR)
LIBS     = -L$(LIBBASE) $(SODIRKEY). $(LIBPATHS) -l$(THISLIB)
SLIBS    = -L$(SLIBBASE) -l$(THISLIB)
SOLABEL  = $(SONAMEKEY)$(SONAME)

OBJS     =  abkassert.o infolines.o paramproc.o abktimer.o platfDepend.o\
	    abkio.o Tausworthe.o r250.o r1279.o seed.o verbosity.o\
	    pathDelims.o abkfunc.o MD5/md5.o crc32.o abkcpunorm.o abkstring.o\
            abkroot.o abkrand_templ.o\
            mmapStreamBuf.o abkmessagebuf.o

SRCS     =  $(OBJS:%.o=%.cxx)

# The four definitions below are for use with the PRINT_VERSION_INFO macro
# defined in ABKCommon/abkversion.h
ABKPACKAGE =`pwd | sed 's/.*\///g'`
ABKSYSTEM  ="`uname -n`\040running\040`uname -sr | sed 's/ /-/g'`"
ABKUSER    =`whoami`
ABKSYMBOLS =-D__ABKPACKAGE__=\"$(ABKPACKAGE)\" -D__ABKSYSTEM__=\"$(ABKSYSTEM)\"\
	-D__ABKUSER__=\"$(ABKUSER)\" -D__ABKLIBBASE__=\"$(LIBBASE)\"

all  : md5_objs lib

md5_objs:
	$(MAKE) -C $(MD5DIR) OPTFLAGS="$(OPTFLAGS)" CCFLAGS="$(CCFLAGS)" all


lib  : md5_objs $(OBJS);
	/bin/rm -f $(ARNAME) $(SODIR)/$(SONAME)
	$(AR) $(ARNAME) $(OBJS)
	$(LD) -o $(SODIR)/$(SONAME) $(OBJS) $(SOLABEL)

TEST_EXES = ABKCommonTest0.exe ABKCommonTest1.exe ABKCommonTest2.exe \
	       ABKCommonTest3.exe ABKCommonTest5.exe ABKCommonTest7.exe \
	       ABKCommonTest9.exe ABKCommonTest10.exe ABKCommonTest11.exe \
	       ABKCommonTest13.exe ABKCommonTest4.exe ABKCommonTest16.exe \
	       ABKCommonTest19.exe ABKCommonTest20.exe ABKCommonTest21.exe \
	       ABKCommonTest6.exe ABKCommonTest22.exe

SANITIZER_FLAGS = -fsanitize=address,undefined -fno-omit-frame-pointer

test : lib $(TEST_EXES)

# Tests link against libABKCommon; without this, "make -jN test" can link a
# test while "lib" is deleting and recreating the library.
$(TEST_EXES): lib


testS: ABKCommonTest0.exeS ABKCommonTest1.exeS ABKCommonTest2.exeS \
       ABKCommonTest3.exeS ABKCommonTest4.exeS ABKCommonTest5.exeS \
       ABKCommonTest6.exeS ABKCommonTest7.exeS ABKCommonTest9.exeS \
       ABKCommonTest10.exeS

.cxx.o:
	$(CC) $(INCLUDES) $(ABKSYMBOLS) -c $*.cxx 

.h.h.gch:
	$(CC) $(INCLUDES) $(ABKSYMBOLS) $*.h

.cxx.exe:
	$(CC) $(INCLUDES) $(ABKSYMBOLS) -c $*.cxx 
	/bin/rm -f $*.exe
	$(LDFINAL) -o $*.exe $*.o $(LIBS)

.cxx.exeS:
	$(CC) $(INCLUDES) $(ABKSYMBOLS) -c $*.cxx 
	/bin/rm -f $*.exeS
	$(LDFINAL) -o $*.exeS $*.o -static $(SLIBS)


numbits:
	/bin/rm -f numbits.exe
	$(CC) -c numbits.cxx
	$(LDFINAL) -o numbits.exe numbits.o

numbitS:
	/bin/rm -f numbits.exeS
	$(CC) -c numbits.cxx
	$(LDFINAL) -static -o numbits.exeS numbits.o

errorcheck:
	/bin/rm -f errorcheck.exe
	$(CC) -c  errorcheck.cxx
	$(LDFINAL) -o errorcheck.exe errorcheck.o

depend:
	/bin/rm -f depend.mk

depend.mk:
	$(CC) $(DEPENDOPT) $(INCLUDES) $(SRCS) >depend.mk

# include depend.mk

clean: md5_clean
	/bin/rm -rf *.o Templates.DB *.exe* TimeTest lib* \
	core cache ir.out MD5/*.o SunWS_cache *notime new.out \
	*.h.gch

md5_clean:
	$(MAKE) -C $(MD5DIR) clean

regression: test
	./regression

sanitizer-regression:
	$(MAKE) clean
	$(MAKE) OPTFLAGS="-g -DABKDEBUG $(SANITIZER_FLAGS)" LDFINAL="g++ $(SANITIZER_FLAGS)" test
	ASAN_OPTIONS=detect_leaks=0 ./ABKCommonTest22.exe /tmp

.PHONY: all lib test testS clean md5_clean regression sanitizer-regression
