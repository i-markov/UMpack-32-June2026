/**************************************************************************
***    
*** Copyright (c) 2026 IonQ, Inc.
*** Copyright (c) 1995-2000 Regents of the University of California,
***               Andrew E. Caldwell, Andrew B. Kahng and Igor L. Markov
*** Copyright (c) 2000-2012 Regents of the University of Michigan,
***               Saurabh N. Adya, Jarrod A. Roy, David A. Papa and
***               Igor L. Markov
***
***  Contact author(s): abk@cs.ucsd.edu, imarkov@umich.edu
***  Original Affiliation:   UCLA, Computer Science Department,
***                          Los Angeles, CA 90095-1596 USA
***
***  Permission is hereby granted, free of charge, to any person obtaining 
***  a copy of this software and associated documentation files (the
***  "Software"), to deal in the Software without restriction, including
***  without limitation 
***  the rights to use, copy, modify, merge, publish, distribute, sublicense, 
***  and/or sell copies of the Software, and to permit persons to whom the 
***  Software is furnished to do so, subject to the following conditions:
***
***  The above copyright notice and this permission notice shall be included
***  in all copies or substantial portions of the Software.
***
*** THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, 
*** EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
*** OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. 
*** IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
*** CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT
*** OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
*** THE USE OR OTHER DEALINGS IN THE SOFTWARE.
***
***
***************************************************************************/


// Created: Mar 18, 2005 by David Papa <iamyou@umich.edu>

#include "ABKCommon/abkcommon.h"
#include "riParams.h"
using std::ostream;
using std::cout;
using std::endl;

RowIroningParameters::RowIroningParameters(int argc, const char *argv[])
 : numPasses(1), windowSize(12), overlap(4), wsChunkCap(3), twoDim(false),
   mixed(false), cluster(false), useDP(true), LRfirst(true),
#ifdef USEFLUTE
   useSteiner(false),
#endif
   smplParams(argc, argv), verb(argc, argv)
{
  NoParams         noParams       (argc, argv);
  BoolParam        helpRequest    ("help", argc, argv);
  BoolParam        helpRequest1   ("h",   argc, argv);
  UnsignedParam numPasses_ ("ironPasses", argc, argv);
  UnsignedParam windowSize_("ironWindow", argc, argv);
  UnsignedParam overlap_   ("ironOverlap", argc, argv);
  UnsignedParam wsChunkCap_("wsChunkCap", argc, argv);
  BoolParam     twoDim_    ("ironTwoDim", argc, argv);
  BoolParam     mixed_     ("ironMixed", argc, argv);
  BoolParam     cluster_   ("ironCluster", argc, argv);
  BoolParam     noDP_      ("ironNoDP", argc, argv);
  BoolParam     RLfirst_   ("ironLeftToRightFirst", argc, argv);
#ifdef USEFLUTE
  BoolParam     useSteiner_("useSteiner", argc, argv);
#endif
  
  if (noParams.found() || helpRequest.found() || helpRequest1.found())
  {
     cout << " Row Ironing Parameters: " 
          << "\n\t -ironPasses <unsigned> "
          << "\n\t -ironWindow <unsigned> "
          << "\n\t -ironOverlap <unsigned> "
	  << "\n\t -ironTwoDim "
	  << "\n\t -ironMixed    (1D and 2D) "
	  << "\n\t -ironCluster "
          << "\n\t -ironNoDP "
	  << "\n\t -ironRightToLeftFirst "
#ifdef USEFLUTE
          << "\n\t -useSteiner "
#endif
	  <<endl;
     return;
  }
  if (numPasses_. found())  numPasses =numPasses_;
  if (windowSize_.found())  windowSize=windowSize_;
  if (overlap_.   found())  overlap   =overlap_;
  if (wsChunkCap_.found())  wsChunkCap=wsChunkCap_;
  if (twoDim_.found()) twoDim = 1;
  if (mixed_.found()) mixed = 1;
  if (cluster_.found()) cluster = 1;
  if (noDP_.found()) useDP = false;
  if (RLfirst_.found()) LRfirst = false;
#ifdef USEFLUTE
  if (useSteiner_.found()) useSteiner = true;
#endif

  smplParams.verb = verb;
  smplParams.verb.setForActions (smplParams.verb.getForActions ()/ 10);
  smplParams.verb.setForMajStats(smplParams.verb.getForMajStats()/ 10);
  smplParams.verb.setForSysRes  (smplParams.verb.getForSysRes  ()/ 10);
}

ostream& operator<<(ostream& out, const RowIroningParameters& params)
{
   out << " Row ironing parameters : "
       << "\n   Num. passes    : " << params.numPasses 
       << "\n   Window size    : " << params.windowSize
       << "\n   Overlap        : " << params.overlap
       << "\n   WS Chunk Cap   : " << params.wsChunkCap
       << "\n   Two Dimensional: " << params.twoDim
       << "\n   Mixed(1D&2D)   : " << params.mixed
       << "\n   Cluster        : " << params.cluster
       << "\n   Use DP         : " << params.useDP
#ifdef USEFLUTE
       << "\n   Use Steiner    : " << params.useSteiner
#endif
       << "\n   Verbosity      : " << params.verb << endl;
   
  return out;
}
