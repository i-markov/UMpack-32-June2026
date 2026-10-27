
# 970821 Created by Max Moroz and Igor Markov,  VLSICAD ABKGroup, UCLA/CS

# 971006 ilm added dependency on Stats
# 990812 ilm upgraded for use with g++  

.SUFFIXES:        .o .cxx .exe .h.gch

# ##################
# standard variables
# ##################

 OPTFLAGS        = -g -DABKDEBUG -DHGDEBUG -DPARTDEBUG -DCAPODEBUG
CCFLAGS          = -fPIC
DEPENDOPT        = -xM1
CC               = g++ $(OPTFLAGS) $(CCFLAGS)
LD               = g++ -shared
LDFINAL          = g++ -Wl,--no-as-needed
AR               = ar cr
SODIRKEY         = -Wl,-rpath,
SONAMEKEY        = -Wl,-soname,
             
UCLA     := $(abspath $(dir $(abspath $(lastword $(MAKEFILE_LIST))))/..)
IMPORT   = NOT_USED
LIBBASE  := $(UCLA)/lib
STL              = 

ARNAME           = lib$(THISLIB).a
SONAME           = lib$(THISLIB).so
SODIR            = $(CURDIR)
SOFULL           = $(SODIR)/$(SONAME)
SOLABEL          = $(SONAMEKEY)$(SONAME)
TESTNAME         = $(THISLIB)Test

# #########################################################
# standard variables that need modification in each project
# #########################################################

THISLIB = ParquetFP
INCLUDES= $(STL) -I$(UCLA)

LIBPATHS         = $(SODIRKEY)$(LIBBASE)
LIBS             =-L. $(SODIRKEY). -l$(THISLIB) $(LIBPATHS) -L$(LIBBASE) -lABKCommon 

BTree_OBJS = basepacking.o btree.o btreecompact.o btreeanneal.o btreeslackeval.o pltobtree.o plsptobtree.o plcompact.o mixedpacking.o baseannealer.o mixedpackingfromdb.o netlist.o

OBJS       = PlToSP.o CommandLine.o FPcommon.o Node.o Nodes.o Net.o Nets.o DB.o ClusterDB.o SeqPair.o SPeval.o AnalytSolve.o Annealer.o SolveMulti.o $(BTree_OBJS)

SRCS       = $(OBJS:%.o=%.cxx)

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

depend:
	/bin/rm -f depend.mk

depend.mk:
	$(CC) $(DEPENDOPT) $(INCLUDES) $(SRCS) *.cxx >depend.mk

clean:        
	rm -rf *.o *.exe *.exe? $(SOFULL) $(ARNAME) core ir.out\
	lib*pure* cache SunWS_cache\
	*notime new.out test.out *.gpl *.png *~ *.h.gch

TEST_EXES = ParquetFPTest0.exe ParquetFPTest1.exe FPEvalTest1.exe FPEvalTest2.exe FPEvalTest3.exe FPEvalTest4.exe
SANITIZER_FLAGS = -fsanitize=address,undefined -fno-omit-frame-pointer

test: lib Parquet.exe $(TEST_EXES)

.PHONY: depend, new, clean, test, lib, regression, sanitizer-regression

# include depend.mk

# #################################
# sample project-specific variables 
# #################################

TEST_EXP        = expected.out
TEST_ARGS       = -f mini2

# ########################
# project-specific targets
# ########################

regression: lib Parquet.exe $(TEST_EXES)
	./regression

sanitizer-regression:
	$(MAKE) clean
	$(MAKE) OPTFLAGS="-g -DABKDEBUG $(SANITIZER_FLAGS)" LDFINAL="g++ $(SANITIZER_FLAGS)" test
	ASAN_OPTIONS=detect_leaks=0 ./ParquetFPTest0.exe TESTS/regression

correct:
	./Parquet.exe $(TEST_ARGS) > $(TEST_EXP)
	perl -n -i -e "print unless /Created|sec/;" $(TEST_EXP)
