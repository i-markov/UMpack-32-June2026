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

// Weighted net-cut evaluator test.
// It assigns explicit edge weights, evaluates the weighted cut, and then
// checks the incremental move update paths.
// Use it to validate weighted partition-cost logic.

#include <iostream>

#include "ABKCommon/paramproc.h"
#include "Partitioning/partitioning.h"
#include "partEvals.h"

using std::cout;
using std::endl;

int main(int argc, char *argv[]) {
  StringParam baseFileName("f", argc, argv);
  abkfatal(baseFileName.found(), " -f baseFileName not found");
  PartitioningProblem problem(baseFileName);

  Partitioning part = problem.getSolnBuffers()[0];

  // NetCutWNetVec eval(problem,part);
  NetCutWNetVec eval(problem);
  eval.resetTo(part);

  cout << eval << endl;

  PartitionIds part0;
  part0.addToPart(0);
  part = Partitioning(part.size(), part0);

  eval.reinitialize();

  cout << " All modules moved to one partition " << endl;
  cout << eval << endl;

  unsigned from = part[4].lowestNumPart();
  PartitionIds old = part[4];
  part[4].setToPart(1);

  //  eval.moveModuleTo(4,from,1);
  eval.moveWithReplication(4, old, part[4]);

  cout << " Module 4 moved to partition 1 " << endl;
  cout << eval << endl;

  from = part[3].lowestNumPart();
  unsigned to = 1 - from;

  cout << "Moving node 3 from " << from << " to " << to << endl;
  part[3].setToPart(to);
  eval.moveModuleTo(3, from, to);
  cout << eval << endl;

  cout << "Moving node back " << endl;
  part[3].setToPart(from);
  eval.moveModuleTo(3, to, from);
  cout << eval << endl;

  return 0;
}
