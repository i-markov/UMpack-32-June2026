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

// Created by Igor Markov
// modified by Andrew Caldwell 08/04/98

// Repeats MLPart over multiple runs and reports best-of-n statistics.

#include "ABKCommon/abkversion.h"
#include "ABKCommon/infolines.h"
#include "ABKCommon/paramproc.h"
#include "PartEvals/partEvals.h"
#include "Partitioning/partProb.h"
#include "Stats/expMins.h"
#include "Stats/trivStats.h"
#include "mlPart.h"

using std::cout;
using std::endl;
using std::vector;

int main(int argc, const char* argv[]) {
  BoolParam helpRequest1("h", argc, argv), helpRequest2("help", argc, argv);
  NoParams noParams(argc, argv);

  StringParam baseFileName("f", argc, argv);
  UnsignedParam numStarts("num", argc, argv);

  UnsignedParam numRuns("runs", argc, argv);
  BoolParam refineSoln("refine", argc, argv);

  PRINT_VERSION_INFO
  cout << CmdLine(argc, argv);

  MLPartParams mlParams(argc, argv);
  cout << endl << mlParams;

  if (helpRequest1.found() || helpRequest2.found() || noParams.found()) {
    cout << " -f   <aux filename>  \n"
         << " -num <multistarts>   \n"
         << " -refine            \n"
         << " -help " << endl;
    return 0;
  }

  abkfatal(baseFileName.found(), " -f baseFileName not found");

  PartitioningProblem::Parameters pParams(argc, argv);
  PartitioningProblem problem(baseFileName, pParams);
  problem.propagateTerminals(mlParams.partFuzziness);
  if (refineSoln.found()) {
    mlParams.seedTopLvlSoln = true;
    problem.setBestSolnNum(0);
  }

  unsigned startsPerRun = 1;
  if (numStarts.found() && !refineSoln.found()) startsPerRun = numStarts;

  problem.reserveBuffers(startsPerRun);

  Partitioning importedInitial(problem.getHGraph().getNumNodes());
  unsigned importedInitialCost = 0;
  UnivPartEval costCheck(PartEvalType::NetCutWNetVec);
  if (refineSoln.found()) {
    importedInitial = problem.getSolnBuffers()[0];
    importedInitialCost = costCheck.computeCost(problem, importedInitial);
    cout << "Initial Net Cut: " << importedInitialCost << endl;
  }

  unsigned totalRuns = 1;
  if (numRuns.found() && !refineSoln.found()) totalRuns = numRuns;

  vector<unsigned> bests;
  vector<double> userTimes;
  vector<double> normUserTimes;
  CPUNormalizer norm;
  double totalRuntime = 0;

  Partitioning bestSeen(problem.getHGraph().getNumNodes());
  unsigned bestSeenCost = UINT_MAX;

  BBPartBitBoardContainer bbBitBoards;
  MaxMem maxMem;

  for (unsigned r = 0; r < totalRuns; r++) {
    cout << "-----   TestRun Number " << r << " of " << totalRuns << " runs \n";

    if (refineSoln.found()) {
      problem.getSolnBuffers()[0] = importedInitial;
      problem.setBestSolnNum(0);
    } else {
      PartitioningBuffer& partBuff = problem.getSolnBuffers();
      unsigned s;
      for (s = partBuff.beginUsedSoln(); s != partBuff.endUsedSoln(); s++) {
        for (unsigned cell = problem.getHGraph().getNumTerminals();
             cell != problem.getHGraph().getNumNodes(); cell++) {
          partBuff[s][cell].setToAll(problem.getNumPartitions());
        }
      }
    }

    problem.propagateTerminals(mlParams.partFuzziness);

    MLPart mlPartitioner(problem, mlParams, bbBitBoards, &maxMem);
    Partitioning& curPart =
        (problem.getSolnBuffers()[problem.getBestSolnNum()]);

    unsigned actualCost = costCheck.computeCost(problem, curPart);
    if (refineSoln.found()) {
      cout << "Final Net Cut: " << actualCost << endl;
      abkfatal(actualCost <= importedInitialCost,
               "Refinement worsened the imported partition");
    }

    if (actualCost < bestSeenCost) {
      bestSeen = curPart;
      bestSeenCost = actualCost;
    }

    cout << "Net Cut: " << actualCost << endl;
    cout << "Part Areas: " << mlPartitioner.getPartitionArea(0) << "  ("
         << mlPartitioner.getPartitionNumNodes(0) << " nodes with terms) "
         << "  /  " << mlPartitioner.getPartitionArea(1) << "  ("
         << mlPartitioner.getPartitionNumNodes(1) << " nodes with terms) "
         << endl;

    bests.push_back(actualCost);
    userTimes.push_back(mlPartitioner.getUserTime());
    totalRuntime += userTimes.back();

    normUserTimes.push_back(mlPartitioner.getUserTime() *
                            norm.getNormalizingFactor());
  }

  cout << "________________ Overall Results ___________________\n";
  cout << "  Best Costs:   " << TrivialStatsWithStdDev(bests) << endl;
  cout << "  RunTime:      " << TrivialStatsWithStdDev(userTimes) << endl;
  cout << "  NRunTime:     " << TrivialStatsWithStdDev(normUserTimes) << endl;

  double runtimePer = totalRuntime / totalRuns;

  if (bests.size() >= 20) {
    ExpectedMins sampler(bests);

    cout << "_______ Average Best Of N/ Runtime _______\n";
    cout << " Ave Best Of 1: " << sampler[1] << "(" << runtimePer << ")"
         << endl;
    cout << " Ave Best Of 2: " << sampler[2] << "(" << runtimePer * 2 << ")"
         << endl;

    if (bests.size() >= 40)
      cout << " Ave Best Of 4: " << sampler[4] << "(" << runtimePer * 4 << ")"
           << endl;

    if (bests.size() >= 50)
      cout << " Ave Best Of 8: " << sampler[8] << "(" << runtimePer * 8 << ")"
           << endl;
  }
  Partitioning& curPart = (problem.getSolnBuffers()[problem.getBestSolnNum()]);
  curPart = bestSeen;
  problem.saveBestSol("best.sol");

  cout << SysInfo();
  cout << " Costs: ";
  for (unsigned i = 0; i < bests.size(); i++) cout << bests[i] << " ";
  cout << endl;
  cout << " Runtime: " << totalRuntime << endl;

  return 0;
}
