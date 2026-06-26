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

// Created by Andrew Caldwell on Apr 27, 1999

// iterative partitioner calls FM several times, with decreasing tolerance
//  if given a problem where there are movable cells if size > 1/2 the
//  original tolerance.

#include "ABKCommon/abkio.h"
#include "iterTolPart.h"

using std::cout;
using std::endl;
using std::max;
using std::ostream;
using std::vector;

IterativeTolPart::IterativeTolPart(PartitioningProblem& problem, MaxMem* maxMem,
                                   Parameters params)
    : _problem(problem),
      _params(params),
      _origMaxCap(problem.getMaxCapacities()),
      _origMinCap(problem.getMinCapacities()) {
  abkfatal(_origMaxCap[0][0] == _origMaxCap[1][0] &&
               _origMinCap[0][0] == _origMinCap[1][0],
           "asymetric max and min values");

  const HGraphFixed& hgraph = _problem.getHGraph();

  const vector<vector<double> >& targets = _problem.getCapacities();
  vector<vector<double> >& maxCaps =
      const_cast<vector<vector<double> >&>(_problem.getMaxCapacities());
  vector<vector<double> >& minCaps =
      const_cast<vector<vector<double> >&>(_problem.getMinCapacities());

  double origTol = _origMaxCap[0][0] - _origMinCap[0][0];

  cout << "Original MaxCap:   " << maxCaps[0][0] << " " << maxCaps[1][0]
       << endl;
  cout << "Original MinCap:   " << minCaps[0][0] << " " << minCaps[1][0] << endl
       << endl;
  cout << "Original Tolerance: " << origTol << endl << endl;

  double maxCellArea = 0;
  itHGFNodeGlobal n;
  for (n = hgraph.nodesBegin(); n != hgraph.nodesEnd(); n++)
    maxCellArea = max(maxCellArea, double(hgraph.getWeight((*n)->getIndex())));

  double ratioTol = _params.initTolRatio * maxCellArea;
  double perctTol =
      (targets[0][0] + targets[1][0]) * _params.initTolPerct / 100.0;

  cout << "Max Cell Area: " << maxCellArea << endl;
  cout << " Ratio:        " << _params.initTolRatio << endl;
  cout << " Ratio Tol:    " << ratioTol << endl;
  cout << " Percent:      " << _params.initTolPerct << endl;
  cout << " Perct Tol:    " << perctTol << endl;

  double newTol = max(ratioTol, max(perctTol, origTol));

  if (newTol <= origTol * 1.01)  // needs to be at least 1% larger, else
  {                              // don't bother
    cout << "____Only Partitioning (w/ orig tol)_______" << endl;
    FMPartitioner fmpart0(_problem, maxMem, _params.fmParams);
    cout << " Cuts: " << _problem.getCosts() << endl;
  } else {
    MultiStartPartitioner::Parameters firstPartParams = _params.fmParams;
    firstPartParams.verb = Verbosity("0 0 0");

    if (_params.firstFMNumPasses != 0)
      firstPartParams.maxNumPasses = _params.firstFMNumPasses;

    if (_params.firstFMUseClip) {
      firstPartParams.useClip = 1;
      if (firstPartParams.maxNumPasses > 1)  // if limited passes,
        firstPartParams.maxNumPasses--;      // do n-1 clip, 1 lifo
    }

    // set the problem to the larget tolerance
    maxCaps[0][0] = targets[0][0] + newTol / 2.0;
    minCaps[0][0] = targets[0][0] - newTol / 2.0;
    maxCaps[1][0] = targets[1][0] + newTol / 2.0;
    minCaps[1][0] = targets[1][0] - newTol / 2.0;
    newTol =
        _problem.getMaxCapacities()[0][0] - _problem.getMinCapacities()[0][0];

    cout << "Resetting tolerance..." << endl;
    cout << " New MaxCap:   " << maxCaps[0][0] << " " << maxCaps[1][0] << endl;
    cout << " New MinCap:   " << minCaps[0][0] << " " << minCaps[1][0] << endl;
    cout << " New Tolerance: " << newTol << endl;

    cout << "____1st Partitioning_______" << endl;
    FMPartitioner fmpart0(_problem, maxMem, firstPartParams);
    cout << " Cuts: " << _problem.getCosts() << endl;

    if (firstPartParams.useClip &&
        _params.firstFMNumPasses > 1) {  // do one lifo pass
      firstPartParams.useClip = 0;
      firstPartParams.maxNumPasses = 1;

      cout << "____Single-Pass LIFO Partitioning_______" << endl;
      FMPartitioner fmpartl(_problem, maxMem, firstPartParams);
      cout << " Cuts: " << _problem.getCosts() << endl;
    }

    // put back the original tolerance
    maxCaps[0][0] = _origMaxCap[0][0];
    maxCaps[1][0] = _origMaxCap[1][0];
    minCaps[0][0] = _origMinCap[0][0];
    minCaps[1][0] = _origMinCap[1][0];
    newTol =
        _problem.getMaxCapacities()[0][0] - _problem.getMinCapacities()[0][0];

    cout << "Orig Tolerance: " << newTol << endl;

    cout << "____2nd Partitioning_______" << endl;
    FMPartitioner fmpart1(_problem, maxMem, params.fmParams);
    cout << " Cuts: " << _problem.getCosts() << endl;
  }
}

ostream& operator<<(ostream& os, const IterativeTolPart::Parameters& params) {
  const char* tfStr[2] = {"false", "true"};

  os << " Iterative Tolerance Partitioner Parameters" << endl;
  os << " initTolRatio:   " << params.initTolRatio << endl;
  os << " firstFMNumPasses: " << params.firstFMNumPasses << endl;
  os << " firstFMUseClip:   " << tfStr[params.firstFMUseClip] << endl;
  os << endl << params.fmParams << endl;

  return os;
}
