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


// FixOr orientation optimizer driver (Bookshelf aux-file path).
//
// Reads a row-based placement from -f <design.aux>, constructs
// OrientOptParameters from any command-line flags (verbosity, pass limit,
// etc.), and runs GreedyOrientOpt to flip/mirror cells so that total
// half-perimeter wirelength decreases.  Reports initial WL, the number
// of cells whose orientation changed, and final WL.  Optionally writes
// the updated placement to -o <out.def>.
//
// This test exercises the pure-RBPlacement code path; it does not require
// the DB or RBPlFromDB packages.

#include <iostream>
#include "ABKCommon/abkcommon.h"
#include "RBPlace/RBPlacement.h"
#include "greedyOrientOpt.h"
using std::cout;
using std::endl;

int main(int argc, const char *argv[])
{
  StringParam auxFileName("f",argc,argv);
  StringParam outFileName("o",argc,argv);
  abkfatal(auxFileName.found(),"Usage: prog -f filename.aux");

  RBPlaceParams rbParam(argc,argv);
  cout << rbParam << endl;
  RBPlacement rbplace(auxFileName, rbParam);

  OrientOptParameters param(argc,argv);
  param.checkSymmetry = false;
  cout << endl << param << endl;

  Timer timer;
  GreedyOrientOpt oo(rbplace,param);
  timer.stop();
  cout << "Greedy Orientation Optimizer took " << timer << endl;
  cout << "Initial WL : " << oo.getInitialWL() << endl
       << "Num Flipped: " << oo.getNumFlips() << endl
       << "Final WL   : " << oo.getFinalWL() << endl;

  if(outFileName.found()) {
        rbplace << oo.getOrients();
        rbplace.savePlacement(outFileName);
  }

  return 0;
}

