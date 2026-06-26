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

#ifndef _ANALYTICALSOLVER_H_
#define _ANALYTICALSOLVER_H_

#include <climits>
#include <map>
#include <unordered_map>

#include <string>
#include <vector>
#include "RBPlace/RBPlacement.h"

class AnalyticalSolver {
 private:
  RBPlacement& _rbplace;
  const PlacementWOrient& _placement;
  std::vector<BBox>&
      _binBBoxs;  // BBox of bin region. used for terminal propagation
  std::vector<std::vector<unsigned> >& _nodes;  // ids of nodes in bin

  // holds final placement and unconstrained solution if CG constraints
  std::vector<std::vector<Point> > _nodesPlacement;
  // holds natural solution to be used only to impose CG constraints
  std::vector<std::vector<Point> > _nodesPlacementNatural;

  std::vector<
      std::unordered_map<unsigned, unsigned, std::hash<unsigned>, std::equal_to<unsigned> > >
      _mapping;
  std::unordered_map<unsigned, unsigned, std::hash<unsigned>, std::equal_to<unsigned> >
      _mappingNodesToBin;

  // sizes of all nodes in the bin. used for imposing CG constraints
  std::vector<double> _binNodesSizes;

  // for gordian L style linearizing in a multi bin scenario
  std::vector<
      std::unordered_map<unsigned, unsigned, std::hash<unsigned>, std::equal_to<unsigned> > >
      _mappingNets;
  std::vector<std::vector<unsigned> > _nets;  // ids of nets in bin

  std::vector<std::vector<double> > _netsLinearizingWeightX;
  std::vector<std::vector<double> > _netsLinearizingWeightY;

  unsigned _numBins;

  unsigned
      _skipNetsLargerThan;  // this is set by default to UINT_MAX
                            // set this with function setSkipNetsLargerThan()

  bool _usePinOffsets;  // set to true by default.
                        // set this with fn setUsePinOffsets()
  double _maxHeightCoreRow;  // if do not use pinoffsets, then pinoffsets are
                             // not used for single row cells

  bool _useLinearizing;  // for gordian style linearizing
  double _avgCellSize;   // for bounding the minimum net length for linearizing

 public:
  AnalyticalSolver();
  AnalyticalSolver(RBPlacement& rbp,
                   std::vector<std::vector<unsigned> >& nodes,
                   std::vector<BBox>& binBBoxs);
  AnalyticalSolver(RBPlacement& rbp, const PlacementWOrient& placement,
                   std::vector<std::vector<unsigned> >& nodes,
                   std::vector<BBox>& binBBoxs);
  ~AnalyticalSolver();

  void setSkipNetsLargerThan(unsigned maxNetDegree);
  void setUsePinOffsets(bool value) { _usePinOffsets = value; }
  void solveQuadraticMin();
  void solveQuadraticMin(double epsilon);
  void combineNaturalUnconstrained();

  void solveSOR(double epsilon, bool naturalSoln = false);
  // do one iteration of SOR in the bin
  double doOneBinPassSOR(unsigned binIdx, bool naturalSoln = false);

  void computeLinearizingWeights();
  void solveLinearMin(unsigned numIter = 5);
  void solveLinearMin(double epsilon, unsigned numIter = 5);
  void solveLinearSOR(double epsilon, bool naturalSoln, unsigned numIter = 5);

  Point getNodeOptLocUnconstrained(unsigned index, unsigned binIdx);
  Point getNodeOptLocNatural(unsigned index, unsigned binIdx);

  void initNodeLocsToBinCenter();
  void initNodeLocs(Point& location, unsigned binIdx);

  std::vector<std::vector<Point> >& getNodeLocs();
  std::vector<std::vector<Point> >& getNodeNaturalLocs();
};

#endif
