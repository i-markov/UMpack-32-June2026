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

#include "ABKCommon/abkassert.h"
#include "ABKCommon/abkversion.h"
#include "ABKCommon/infolines.h"
#include "ABKCommon/paramproc.h"
#include "PartEvals/partEvals.h"
#include "fmPart.h"

#include <iostream>

using std::cout;
using std::endl;

// Basic FM partitioner smoke test.
// Loads the benchmark named by -f, builds an FMPartitioner, and reports
// nonzero-net statistics for the resulting best solution.
// With -refine, it starts from the imported .sol solution and checks that
// FM does not worsen it.
// Optionally writes the best partition to best.sol when -save is given.
int main(int argc, const char* argv[]) {
  BoolParam helpRequest1("h", argc, argv);
  BoolParam helpRequest2("help", argc, argv);

  StringParam baseFileName("f", argc, argv);
  UnsignedParam numStarts("num", argc, argv);

  BoolParam saveSoln("save", argc, argv);
  BoolParam refineSoln("refine", argc, argv);

  PRINT_VERSION_INFO

  PartitioningProblem::Parameters params(argc, argv);
  PartitionerParams pParams(argc, argv);
  if (helpRequest1.found() || helpRequest2.found()) return 0;
  abkfatal(baseFileName.found(), " -f baseFileName not found");
  cout << params << endl;
  PartitioningProblem problem(baseFileName, params);
  if (!saveSoln.found()) problem.discardHGraphNames();
  problem.propagateTerminals();
  if (refineSoln.found()) problem.setBestSolnNum(0);
  if (numStarts.found() && !refineSoln.found())
    problem.reserveBuffers(numStarts);
  cout << pParams << endl;
  MaxMem maxMem;
  UniversalPartEval eval(PartEvalType::NetCut2wayWWeights);
  if (refineSoln.found()) {
    const Partitioning& initialPart = problem.getBestSoln();
    unsigned initialCost = eval.computeCost(problem, initialPart);
    cout << "Initial Net Cut: " << initialCost << endl;

    FMPartitioner fm(problem, &maxMem, pParams, true, true);
    fm.runMultiStart();

    const Partitioning& finalPart = problem.getBestSoln();
    unsigned finalCost = eval.computeCost(problem, finalPart);
    cout << "Final Net Cut: " << finalCost << endl;
    abkfatal(finalCost <= initialCost,
             "Refinement worsened the imported partition");
  } else {
    FMPartitioner fm(problem, &maxMem, pParams);
  }
  eval.printStatsForNon0CostNets(problem);

  if (saveSoln.found()) problem.saveBestSol("best.sol");

  problem.printLargestCellStats();

  cout << SysInfo();

  return 0;
}
