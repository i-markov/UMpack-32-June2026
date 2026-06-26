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

#ifndef THUNDERPART_H
#define THUNDERPART_H

#include "ABKCommon/abkrand.h"
#include "Geoms/bbox.h"
#include "Geoms/point.h"
#include "HGraphWDims/hgWDims.h"
#include "Placement/placeWOri.h"
#include "RBPlace/RBPlacement.h"
#include "capoBin.h"

/*
 * <aaronnn> ThunderPart does global placement refinement
 *           by relocating cells to different bins; optimizes HPWL
   WARNING:  This code currently improves intermediate wirelength
   but does not improve final wirelength.  It is not clear why,
   so use this code at your own risk.
 */

class ThunderPartParams {
 public:
  Verbosity verb;

  ThunderPartParams() {}
  ThunderPartParams(const ThunderPartParams& srcParams) {
    verb = srcParams.verb;
  }
};

class ThunderPart {
 public:
  ThunderPart(std::vector<CapoBin*>& bins, const HGraphWDimensions& hgraph,
              const PlacementWOrient& placement, const RBPlacement& rbplace,
              ThunderPartParams params);
  void partitionize(unsigned);
  double getHPWL();
  const std::vector<std::vector<unsigned> >& getResultBins() {
    return _resultBins;
  }

 private:
  /*
   * private functions
   */
  void initialize();
  void recomputeExposedNets();
  void updateExposedNet(unsigned edgeIdx);
  void updateHPWL(unsigned edgeID);
  bool netIsCut(unsigned) const;
  void moveCell(unsigned nodeID, unsigned binID, Point newLoc);
  Point findBestCaseDest(unsigned vertexIdx, unsigned binMovedTo);
  bool noVisibleEdges();
  const HGFEdge& selectVisibleEdge();
  unsigned getWeberBin(unsigned nodeIdx);
  bool willOverfillBin(unsigned nodeID, unsigned binID);
  bool isOutOfCore(Point& point);
  bool isMovable(unsigned nodeID);
  void announceHPWL(char c, double HPWL);

  /*
   * algorizms
   */
  void simpleGreedy(unsigned numMoves);
  // void hillClimb(unsigned numMoves);

  /*
   * input
   */
  const std::vector<CapoBin*>& _bins;
  const RBPlacement& _rbplace;
  const HGraphWDimensions& _hgraph;

  /*
   * output
   */
  std::vector<std::vector<unsigned> > _resultBins;

  /*
   * runnning (incrementally updated) values
   */
  PlacementWOrient _placement;  // location of cells, placed by ThunderPart
  std::vector<unsigned>
      _nodeBins;               // cellID -> binID, NO_BIN to indicate no bin
  std::vector<double> _HPWL;  // used to track HPWL for each net
  std::vector<double> _binCellAreas;  // total area of cells in bin

  // double _sumNetWeights;
  double _totalHPWL;
  ThunderPartParams _params;

  /*
   * stuff for edge selection
   */
  std::vector<std::vector<unsigned> >
      _edgeSelectionBuckets;  // put edges in buckets, by their edge degrees
  std::vector<unsigned>
      _edgeSelectionBucketLoc;  // idx of edge in an edge-selection-bucket

  /*
   * misc
   */
  RandomRawUnsigned _randomNum;
  bool coinFlip();
  bool coinFlip(double sideProb);
};

#endif /* THUNDERPART_H */
