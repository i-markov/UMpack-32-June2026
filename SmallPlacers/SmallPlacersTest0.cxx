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

// SmallPlacers branch-and-bound regression driver.
// It runs the small placement solver multiple times on one .aux design to
// exercise the bound computation, repeated search, and optional plotting.
// Use -f <auxFile>, -num <runs>, -plot <file>, and -save <file> to control it.

// CHANGES

#include <cstdlib>
#include <iostream>

#include "ABKCommon/paramproc.h"
#include "ABKCommon/abktimer.h"
#include "ABKCommon/infolines.h"
#include "SmallPlacement/smallPlaceProb.h"
#include "sRowBBPlacer.h"
#include "sRowSmPlacer.h"

using std::cout;
using std::endl;

int main(int argc, const char* argv[]) {
  StringParam auxFile("f", argc, argv);
  UnsignedParam num("num", argc, argv);
  StringParam plotFile("plot", argc, argv);
  StringParam save("save", argc, argv);
  BoolParam helpRequest0("h", argc, argv);
  BoolParam helpRequest1("help", argc, argv);

  SmallPlParams smParams(argc, argv);
  if (helpRequest0.found() || helpRequest1.found()) exit(-1);

  cout << smParams;

  abkfatal(auxFile.found(), "-f <auxFile> not found");

  SmallPlacementBitBoardContainer bitBoards;

  SmallPlacementProblem problem(auxFile, bitBoards);

  unsigned numTests = (num.found()) ? num : 10;

  MaxMem maxMem;

  Timer t;

  for (unsigned i = 0; i < numTests; i++) {
    // BranchSmallPlacer placer(problem, smParams);
    SmallPlacer placer(problem, smParams, &maxMem);
    problem.getSolnBuffer()[0].populated = false;
    smParams.verb = Verbosity("0_0_0");
  }

  t.stop();

  if (plotFile.found()) problem.plot(plotFile);

  Placement& placement = problem.getBestSoln().placement;
  cout << "Final WL is " << problem.calculateWL(placement) << endl;
  cout << numTests << " placements took " << t.getUserTime() << "  "
       << t.getUserTime() / (double)numTests << " each" << endl;

  if (save.found()) problem.save(save);

  return 0;
}
