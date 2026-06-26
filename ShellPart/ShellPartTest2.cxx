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

// Iterative tolerance partitioning driver.
// It loads a partitioning problem from -f, reserves buffers for multiple
// starts, runs the iterative tolerance partitioner, and prints runtime and
// expected-minimum summaries.

#include "PartEvals/partEvals.h"
#include "ABKCommon/abkcpunorm.h"
#include "ABKCommon/infolines.h"
#include "Partitioning/partProb.h"
#include "Stats/expMins.h"
#include "Stats/trivStats.h"
#include "iterTolPart.h"

using std::cout;
using std::endl;

int main(int argc, const char *argv[]) {
  StringParam baseFileName("f", argc, argv);
  UnsignedParam num("num", argc, argv);

  PartitioningProblem::Parameters params(argc, argv);

  abkfatal(baseFileName.found(), " -f baseFileName not found");
  PartitioningProblem problem(baseFileName, params);

  unsigned numStarts = 1;
  if (num.found()) numStarts = num;

  problem.reserveBuffers(numStarts);

  problem.propagateTerminals();

  IterativeTolPart::Parameters iTolPartParams(argc, argv);
  cout << iTolPartParams << endl << endl;

  MaxMem maxMem;
  Timer iterTime;
  IterativeTolPart partitioner(problem, &maxMem, iTolPartParams);
  iterTime.stop();

  double runtime = iterTime.getUserTime();
  CPUNormalizer norm;
  double nruntime = runtime * norm.getNormalizingFactor();

  cout << "_______________Final Results_________________" << endl;
  cout << " Cuts: " << TrivialStats(problem.getCosts()) << endl;
  cout << " RunTime:  " << runtime / numStarts << " sec per start" << endl;
  cout << " NRunTime: " << runtime / numStarts << " sec per start" << endl;

  ExpectedMins sampler(problem.getCosts());
  cout << " Ave Best of 1: " << sampler[1] << "(" << nruntime / numStarts << ")"
       << endl;
  cout << " Ave Best of 2: " << sampler[2] << "(" << 2 * (nruntime / numStarts)
       << ")" << endl;
  cout << " Ave Best of 4: " << sampler[4] << "(" << 4 * (nruntime / numStarts)
       << ")" << endl;
  cout << " Ave Best of 8: " << sampler[8] << "(" << 8 * (nruntime / numStarts)
       << ")" << endl;

  cout << SysInfo();
  return 0;
}
