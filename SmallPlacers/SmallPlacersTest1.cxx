/**************************************************************************
***
*** Copyright (c) 2026 IonQ, Inc.
*** Copyright (c) 1995-2000 Regents of the University of California,
***               Andrew E. Caldwell, Andrew B. Kahng and Igor L. Markov
*** Copyright (c) 2000-2012 Regents of the University of Michigan,
***               Saurabh N. Adya, Jarrod A. Roy, David A. Papa and
***               Igor L. Markov
***
***  Contact author(s): abk@cs.ucsd.edu, igor.markov1@gmail.com
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

// SmallPlacers timing and bound comparison driver.
// It first computes a branch-and-bound bound, then repeats the search
// several times to compare runtime and final wirelength on the same design.
// Use -f <auxFile> and -num <runs> to tune the workload.

// CHANGES

#include <iostream>

#include "ABKCommon/abktimer.h"
#include "SmallPlacement/smallPlaceProb.h"
#include "sRowBBPlacer.h"

int main(int argc, const char *argv[]) {
  StringParam auxFile("f", argc, argv);
  UnsignedParam num("num", argc, argv);

  SmallPlParams smParams(argc, argv);

  abkfatal(auxFile.found(), "-f <auxFile> not found");

  SmallPlacementProblem problem(auxFile);

  unsigned numTests = (num.found()) ? num : 10;

  cout << " First run to compute the bound " << endl;
  BranchSmallPlacer pl1(problem, smParams);
  cout << " First run done " << endl;

  cout << " Running " << numTests << " runs "
       << " (at most one with high verbosity) " << endl;

  Timer t;

  for (unsigned i = 0; i < numTests; i++) {
    BranchSmallPlacer placer(problem, smParams);
    //	problem.getSolnBuffer()[0].populated = false;
    smParams.verb = Verbosity("0_0_0");
  }

  t.stop();

  BaseSmallPlacer placer(problem, smParams);
  cout << "Final WL is " << placer.calculateWL() << endl;
  cout << numTests << " placements took " << t.getUserTime() << "  "
       << t.getUserTime() / (double)numTests << " each" << endl;

  return 0;
}
