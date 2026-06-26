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

// baseSmallPl.h
// Base Small Placer
//    created by Andrew Caldwell  981206

#ifndef __SMALLPL_BASE_H__
#define __SMALLPL_BASE_H__

#include "ABKCommon/infolines.h"
#include "SmallPlacement/smallPlaceProb.h"
// #include "spNetStacks.h"

class SmallPlacerNetStacks;

class BaseSmallPlacer {
 public:
  struct Parameters {
    enum SmPlAlgo { BranchAndBound, DynamicProgram };
    // bnb enums.
    enum MinWL { None, HalfRowHeight, RowHeight, TwiceRowHeight };
    enum AssignMethod { LeftOnly, LeftRight };

    // dp enums.
    enum StartMethod { RandomStart, SmartStart };
    enum SplitMethod { EvenSplit, RandomSplit, AllSplit };

    // general parameters.
    Verbosity verb;
    bool useInitSoln;  // use the init soln, if it's given one
    SmPlAlgo algo;     // which algo to use.

    // bnb parameters.
    AssignMethod assignMethod;
    MinWL minWL;
    bool subsetBound;

    // dp parameters.
    StartMethod startSoln;
    SplitMethod splitType;

    Parameters()
        : verb("0_0_0"),
          useInitSoln(true),
          algo(BranchAndBound),
          assignMethod(LeftRight),
          minWL(None),
          subsetBound(false),
          startSoln(SmartStart),
          splitType(AllSplit) {}

    Parameters(int argc, const char* argv[]);

    friend std::ostream& operator<<(std::ostream& os, const Parameters& params);
  };

 protected:
  SmallPlacementProblem& _problem;

  //    SmallPlacerNetStacks*   _netStacks;

  Parameters _params;

  MaxMem* _maxMem;

  // just a handy function..too slow to be used during optimization..
  double calculateWL(const Placement& placement);

 public:
  BaseSmallPlacer(SmallPlacementProblem& problem, Parameters params,
                  MaxMem* maxMem)
      : _problem(problem), _params(params), _maxMem(maxMem) {}
  virtual ~BaseSmallPlacer() {}

  double calculateWL() { return calculateWL(_problem.getBestSoln().placement); }
};

typedef BaseSmallPlacer::Parameters SmallPlParams;

#endif
