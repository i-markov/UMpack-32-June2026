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

#define SOR_FACTOR 1.7

#include "analytPl.h"

#include <algorithm>

using std::cout;
using std::endl;
using std::max;
using std::vector;

AnalyticalSolver::~AnalyticalSolver() {
  _mapping.clear();
  _mappingNodesToBin.clear();
  _nodesPlacement.clear();
  _nodesPlacementNatural.clear();
  _binNodesSizes.clear();
  _mappingNets.clear();
  _nets.clear();
  _netsLinearizingWeightX.clear();
  _netsLinearizingWeightY.clear();
}

AnalyticalSolver::AnalyticalSolver(RBPlacement& rbp,
                                   vector<vector<unsigned> >& nodes,
                                   vector<BBox>& binBBoxs)
    : _rbplace(rbp),
      _placement(rbp.getPlacement()),
      _binBBoxs(binBBoxs),
      _nodes(nodes),
      _skipNetsLargerThan(UINT_MAX),
      _usePinOffsets(true),
      _useLinearizing(false),
      _avgCellSize(0.0) {
  _numBins = _nodes.size();
  _nodesPlacement.resize(_numBins);
  _nodesPlacementNatural.resize(_numBins);
  _mapping.resize(_numBins);
  _mappingNets.resize(_numBins);
  _nets.resize(_numBins);
  _netsLinearizingWeightX.resize(_numBins);
  _netsLinearizingWeightY.resize(_numBins);

  _binNodesSizes.resize(_numBins);

  const HGraphWDimensions& hgWDims = _rbplace.getHGraph();

  double binNodesSize = 0;
  _avgCellSize = 0.0;
  unsigned numCells = 0;

  for (unsigned i = 0; i < _numBins; ++i) {
    unsigned numEdgesInBin = 0;

    binNodesSize = 0;
    _nodesPlacement[i].resize(_nodes[i].size());
    _nodesPlacementNatural[i].resize(_nodes[i].size());
    for (unsigned j = 0; j < _nodes[i].size(); ++j) {
      unsigned index = _nodes[i][j];
      _mapping[i][index] = j;
      _nodesPlacement[i][j].x = _placement[index].x;
      _nodesPlacement[i][j].y = _placement[index].y;
      _nodesPlacementNatural[i][j] = _nodesPlacement[i][j];
      _mappingNodesToBin[index] = i;

      binNodesSize +=
          hgWDims.getNodeWidth(index) * hgWDims.getNodeHeight(index);
      ++numCells;

      // get the net info now
      const HGFNode& origNode = hgWDims.getNodeByIdx(index);
      itHGFEdgeLocal e;
      for (e = origNode.edgesBegin(); e != origNode.edgesEnd(); ++e) {
        HGFEdge& edge = **e;
        unsigned edgeId = edge.getIndex();
        _mappingNets[i][edgeId] = numEdgesInBin;
        _nets[i].push_back(edgeId);
        ++numEdgesInBin;
      }
    }
    _binNodesSizes[i] = binNodesSize;
    _avgCellSize += binNodesSize;

    _netsLinearizingWeightX[i].resize(numEdgesInBin, 1.0);
    _netsLinearizingWeightY[i].resize(numEdgesInBin, 1.0);
  }
  _avgCellSize /= numCells;
  _avgCellSize = 2 * sqrt(_avgCellSize);
  _maxHeightCoreRow = _rbplace.getMaxHeightCoreRow();
}

AnalyticalSolver::AnalyticalSolver(RBPlacement& rbp,
                                   const PlacementWOrient& placement,
                                   vector<vector<unsigned> >& nodes,
                                   vector<BBox>& binBBoxs)
    : _rbplace(rbp),
      _placement(placement),
      _binBBoxs(binBBoxs),
      _nodes(nodes),
      _skipNetsLargerThan(UINT_MAX),
      _usePinOffsets(true),
      _useLinearizing(false),
      _avgCellSize(0.0) {
  _numBins = _nodes.size();
  _nodesPlacement.resize(_numBins);
  _nodesPlacementNatural.resize(_numBins);
  _mapping.resize(_numBins);
  _mappingNets.resize(_numBins);
  _nets.resize(_numBins);
  _netsLinearizingWeightX.resize(_numBins);
  _netsLinearizingWeightY.resize(_numBins);

  _binNodesSizes.resize(_numBins);

  const HGraphWDimensions& hgWDims = _rbplace.getHGraph();

  double binNodesSize = 0;
  _avgCellSize = 0.0;
  unsigned numCells = 0;

  for (unsigned i = 0; i < _numBins; ++i) {
    unsigned numEdgesInBin = 0;

    binNodesSize = 0;
    _nodesPlacement[i].resize(_nodes[i].size());
    _nodesPlacementNatural[i].resize(_nodes[i].size());
    for (unsigned j = 0; j < _nodes[i].size(); ++j) {
      unsigned index = _nodes[i][j];
      _mapping[i][index] = j;
      _nodesPlacement[i][j].x = _placement[index].x;
      _nodesPlacement[i][j].y = _placement[index].y;
      _nodesPlacementNatural[i][j] = _nodesPlacement[i][j];
      _mappingNodesToBin[index] = i;

      binNodesSize +=
          hgWDims.getNodeWidth(index) * hgWDims.getNodeHeight(index);
      ++numCells;

      // get the net info now
      const HGFNode& origNode = hgWDims.getNodeByIdx(index);
      itHGFEdgeLocal e;
      for (e = origNode.edgesBegin(); e != origNode.edgesEnd(); ++e) {
        HGFEdge& edge = **e;
        unsigned edgeId = edge.getIndex();
        _mappingNets[i][edgeId] = numEdgesInBin;
        _nets[i].push_back(edgeId);
        ++numEdgesInBin;
      }
    }
    _binNodesSizes[i] = binNodesSize;
    _avgCellSize += binNodesSize;

    _netsLinearizingWeightX[i].resize(numEdgesInBin, 1.0);
    _netsLinearizingWeightY[i].resize(numEdgesInBin, 1.0);
  }
  _avgCellSize /= numCells;
  _avgCellSize = 2 * sqrt(_avgCellSize);
  _maxHeightCoreRow = _rbplace.getMaxHeightCoreRow();
}

void AnalyticalSolver::setSkipNetsLargerThan(unsigned maxNetDegree) {
  _skipNetsLargerThan = maxNetDegree;
}

vector<vector<Point> >& AnalyticalSolver::getNodeLocs() {
  return _nodesPlacement;
}

vector<vector<Point> >& AnalyticalSolver::getNodeNaturalLocs() {
  return _nodesPlacementNatural;
}

void AnalyticalSolver::initNodeLocsToBinCenter() {
  for (unsigned b = 0; b < _numBins; ++b) {
    Point binCenter = _binBBoxs[b].getGeomCenter();
    initNodeLocs(binCenter, b);
  }
}

void AnalyticalSolver::initNodeLocs(Point& location, unsigned binIdx) {
  for (unsigned i = 0; i < _nodes[binIdx].size(); ++i) {
    _nodesPlacement[binIdx][i].x = location.x;
    _nodesPlacement[binIdx][i].y = location.y;
    _nodesPlacementNatural[binIdx][i].x = location.x;
    _nodesPlacementNatural[binIdx][i].y = location.y;
  }
}

void AnalyticalSolver::solveLinearMin(unsigned numIter) {
  _useLinearizing = true;
  for (unsigned i = 0; i < numIter; ++i) {
    solveQuadraticMin();
    computeLinearizingWeights();
  }
}

void AnalyticalSolver::solveLinearMin(double epsilon, unsigned numIter) {
  _useLinearizing = true;
  for (unsigned i = 0; i < numIter; ++i) {
    solveQuadraticMin(epsilon);
    computeLinearizingWeights();
  }
}

void AnalyticalSolver::solveLinearSOR(double epsilon, bool naturalSoln,
                                      unsigned numIter) {
  _useLinearizing = true;
  for (unsigned i = 0; i < numIter; ++i) {
    solveSOR(epsilon, naturalSoln);
    computeLinearizingWeights();
  }
}

void AnalyticalSolver::solveQuadraticMin() {
  double maxDimension = max(_binBBoxs[0].getHeight(), _binBBoxs[0].getWidth());
  double epsilon = (maxDimension * sqrt(double(_numBins))) / 100;
  solveQuadraticMin(epsilon);
}

void AnalyticalSolver::solveQuadraticMin(double epsilon) {
  if (_rbplace.getParameters().verb.getForMajStats() > 1)
    cout << endl << "SOR Solver - Getting Unconstrained Solution ..." << endl;
  solveSOR(epsilon, 0);  // get unconstrained solution
  if (_rbplace.getParameters().verb.getForMajStats() > 1)
    cout << endl << "SOR Solver - Getting Natural Solution ..." << endl;
  solveSOR(epsilon, 1);  // get natural solution
  combineNaturalUnconstrained();
}

void AnalyticalSolver::combineNaturalUnconstrained() {
  const HGraphWDimensions& hgWDims = _rbplace.getHGraph();
  for (unsigned b = 0; b < _numBins; ++b) {
    double xLambda;
    double yLambda;
    double xWeightUnconstrained = 0;
    double xWeightNatural = 0;
    double yWeightUnconstrained = 0;
    double yWeightNatural = 0;
    for (unsigned i = 0; i < _nodes[b].size(); ++i) {
      unsigned index = _nodes[b][i];
      double nodeSize =
          hgWDims.getNodeWidth(index) * hgWDims.getNodeHeight(index);
      xWeightUnconstrained += nodeSize * _nodesPlacement[b][i].x;
      yWeightUnconstrained += nodeSize * _nodesPlacement[b][i].y;
      xWeightNatural += nodeSize * _nodesPlacementNatural[b][i].x;
      yWeightNatural += nodeSize * _nodesPlacementNatural[b][i].y;
    }
    Point CG = _binBBoxs[b].getGeomCenter();
    xLambda =
        (_binNodesSizes[b] * CG.x - xWeightUnconstrained) / xWeightNatural;
    yLambda =
        (_binNodesSizes[b] * CG.y - yWeightUnconstrained) / yWeightNatural;

    for (unsigned i = 0; i < _nodes[b].size(); ++i) {
      _nodesPlacement[b][i].x += xLambda * _nodesPlacementNatural[b][i].x;
      _nodesPlacement[b][i].y += yLambda * _nodesPlacementNatural[b][i].y;
    }
  }
}

void AnalyticalSolver::solveSOR(double epsilon, bool naturalSoln) {
  double change = DBL_MAX;
  unsigned numIter = 0;

  while (change > epsilon && numIter < 500) {
    change = 0;
    numIter++;
    for (unsigned i = 0; i < _numBins; ++i)
      change = max(doOneBinPassSOR(i, naturalSoln), change);
    if (_rbplace.getParameters().verb.getForMajStats() > 1)
      cout << "\tIteration " << numIter << "\t Max Movement " << change << endl;
  }
  if (_rbplace.getParameters().verb.getForActions() > 1)
    cout << "Epsilon is " << epsilon << "  numIter  " << numIter << endl;
}

double AnalyticalSolver::doOneBinPassSOR(unsigned binIdx, bool naturalSoln) {
  double xchange, ychange, xOvershoot, yOvershoot, indChange;
  double change = 0;
  for (unsigned i = 0; i < _nodes[binIdx].size(); ++i) {
    unsigned index = _nodes[binIdx][i];
    Point newLoc;
    if (naturalSoln == 0)  // unconstrained solution
    {
      newLoc = getNodeOptLocUnconstrained(index, binIdx);
      xchange = newLoc.x - _nodesPlacement[binIdx][i].x;
      ychange = newLoc.y - _nodesPlacement[binIdx][i].y;
      xOvershoot = xchange * SOR_FACTOR;
      yOvershoot = ychange * SOR_FACTOR;
      _nodesPlacement[binIdx][i].x += xOvershoot;
      _nodesPlacement[binIdx][i].y += yOvershoot;
    } else  // natural solution
    {
      newLoc = getNodeOptLocNatural(index, binIdx);
      xchange = newLoc.x - _nodesPlacementNatural[binIdx][i].x;
      ychange = newLoc.y - _nodesPlacementNatural[binIdx][i].y;
      xOvershoot = xchange * SOR_FACTOR;
      yOvershoot = ychange * SOR_FACTOR;
      _nodesPlacementNatural[binIdx][i].x += xOvershoot;
      _nodesPlacementNatural[binIdx][i].y += yOvershoot;
    }
    indChange = fabs(xchange) + fabs(ychange);
    change = max(change, indChange);
  }
  return change;
}

Point AnalyticalSolver::getNodeOptLocUnconstrained(unsigned index,
                                                   unsigned binIdx) {
  Point optLoc;
  const HGraphWDimensions& hgWDims = _rbplace.getHGraph();
  const HGFNode& origNode = hgWDims.getNodeByIdx(index);
  itHGFEdgeLocal e;
  itHGFNodeLocal v;
  double xSum = 0, ySum = 0;
  unsigned edgeOffset = 0;
  double denominatorX = 0;
  double denominatorY = 0;
  for (e = origNode.edgesBegin(); e != origNode.edgesEnd(); ++e, ++edgeOffset) {
    HGFEdge& edge = **e;
    unsigned edgeId = edge.getIndex();

    Point pOffsetOrigNode(0, 0);
    if (_usePinOffsets || hgWDims.getNodeHeight(index) > _maxHeightCoreRow)
      pOffsetOrigNode = hgWDims.edgesOnNodePinOffset(edgeOffset, index,
                                                     _rbplace.getOrient(index));

    unsigned edgeDegree = edge.getDegree();
    double netXSum = 0;
    double netYSum = 0;
    unsigned nodeOffset = 0;

    Point pOffset(0, 0);
    if (edgeDegree <= _skipNetsLargerThan) {
      for (v = (*e)->nodesBegin(); v != (*e)->nodesEnd(); v++, nodeOffset++) {
        HGFNode& node = (**v);
        unsigned nodeIndex = node.getIndex();

        if (_usePinOffsets ||
            hgWDims.getNodeHeight(nodeIndex) > _maxHeightCoreRow)
          pOffset = hgWDims.nodesOnEdgePinOffset(nodeOffset, edgeId,
                                                 _rbplace.getOrient(nodeIndex));
        if (_mapping[binIdx].count(nodeIndex) == 0)  // node not in this
                                                     // bin. so treat as fixed
        {
          double pinXLoc, pinYLoc;
          // node not in all bins
          if (_mappingNodesToBin.count(nodeIndex) == 0) {
            pinXLoc = _placement[nodeIndex].x + pOffset.x - pOffsetOrigNode.x;
            pinYLoc = _placement[nodeIndex].y + pOffset.y - pOffsetOrigNode.y;
          } else {
            unsigned nodeBinId = _mappingNodesToBin[nodeIndex];
            unsigned mappedIdx = _mapping[nodeBinId][nodeIndex];
            pinXLoc = _nodesPlacement[nodeBinId][mappedIdx].x + pOffset.x -
                      pOffsetOrigNode.x;
            pinYLoc = _nodesPlacement[nodeBinId][mappedIdx].y + pOffset.y -
                      pOffsetOrigNode.y;
          }
          // do terminal propagation
          if (pinXLoc < _binBBoxs[binIdx].xMin)
            pinXLoc = _binBBoxs[binIdx].xMin;
          else if (pinXLoc > _binBBoxs[binIdx].xMax)
            pinXLoc = _binBBoxs[binIdx].xMax;

          if (pinYLoc < _binBBoxs[binIdx].yMin)
            pinYLoc = _binBBoxs[binIdx].yMin;
          else if (pinYLoc > _binBBoxs[binIdx].yMax)
            pinYLoc = _binBBoxs[binIdx].yMax;

          netXSum += pinXLoc;
          netYSum += pinYLoc;
        } else if (nodeIndex != index) {
          unsigned mappedIdx = _mapping[binIdx][nodeIndex];
          netXSum += _nodesPlacement[binIdx][mappedIdx].x + pOffset.x -
                     pOffsetOrigNode.x;
          netYSum += _nodesPlacement[binIdx][mappedIdx].y + pOffset.y -
                     pOffsetOrigNode.y;
        }
      }
    }
    if (edgeDegree > 0 && edgeDegree <= _skipNetsLargerThan) {
      double netWt = 2.0 / edgeDegree;
      if (_useLinearizing) {
        unsigned binEdgeId = _mappingNets[binIdx][edgeId];
        netXSum = netWt * netXSum * _netsLinearizingWeightX[binIdx][binEdgeId];
        netYSum = netWt * netYSum * _netsLinearizingWeightY[binIdx][binEdgeId];
        denominatorX += (netWt * (edgeDegree - 1.0) *
                         _netsLinearizingWeightX[binIdx][binEdgeId]);
        denominatorY += (netWt * (edgeDegree - 1.0) *
                         _netsLinearizingWeightY[binIdx][binEdgeId]);
      } else {
        netXSum = netWt * netXSum;
        netYSum = netWt * netYSum;
        denominatorX += (netWt * (edgeDegree - 1.0));
        denominatorY += (netWt * (edgeDegree - 1.0));
      }
      xSum += netXSum;
      ySum += netYSum;
    }
  }

  if (denominatorX != 0 && denominatorY != 0) {
    xSum /= denominatorX;
    ySum /= denominatorY;
    optLoc.x = xSum;
    optLoc.y = ySum;
  } else
    optLoc = _nodesPlacement[binIdx][_mapping[binIdx][index]];

  return optLoc;
}

Point AnalyticalSolver::getNodeOptLocNatural(unsigned index, unsigned binIdx) {
  Point optLoc;
  const HGraphWDimensions& hgWDims = _rbplace.getHGraph();
  const HGFNode& origNode = hgWDims.getNodeByIdx(index);
  itHGFEdgeLocal e;
  itHGFNodeLocal v;
  double xSum = 0, ySum = 0;
  unsigned edgeOffset = 0;
  double denominatorX = 0;
  double denominatorY = 0;
  for (e = origNode.edgesBegin(); e != origNode.edgesEnd(); ++e, ++edgeOffset) {
    HGFEdge& edge = **e;
    unsigned edgeId = edge.getIndex();

    Point pOffsetOrigNode(0, 0);
    if (_usePinOffsets)
      pOffsetOrigNode = hgWDims.edgesOnNodePinOffset(edgeOffset, index,
                                                     _rbplace.getOrient(index));

    unsigned edgeDegree = edge.getDegree();
    double netXSum = 0;
    double netYSum = 0;
    unsigned nodeOffset = 0;
    Point pOffset(0, 0);

    if (edgeDegree <= _skipNetsLargerThan) {
      for (v = (*e)->nodesBegin(); v != (*e)->nodesEnd(); v++, nodeOffset++) {
        HGFNode& node = (**v);
        unsigned nodeIndex = node.getIndex();

        if (_usePinOffsets)
          pOffset = hgWDims.nodesOnEdgePinOffset(nodeOffset, edgeId,
                                                 _rbplace.getOrient(nodeIndex));

        if (_mapping[binIdx].count(nodeIndex) == 0)  // node not in this
                                                     // bin. so treat as fixed
        {
          // treat all fixed objects fixed to 0
        } else if (nodeIndex != index) {
          unsigned mappedIdx = _mapping[binIdx][nodeIndex];
          netXSum += _nodesPlacement[binIdx][mappedIdx].x + pOffset.x -
                     pOffsetOrigNode.x;
          netYSum += _nodesPlacement[binIdx][mappedIdx].y + pOffset.y -
                     pOffsetOrigNode.y;
        }
      }
    }
    if (edgeDegree > 0 && edgeDegree <= _skipNetsLargerThan) {
      double netWt = 2.0 / edgeDegree;
      if (_useLinearizing) {
        unsigned binEdgeId = _mappingNets[binIdx][edgeId];
        netXSum = netWt * netXSum * _netsLinearizingWeightX[binIdx][binEdgeId];
        netYSum = netWt * netYSum * _netsLinearizingWeightY[binIdx][binEdgeId];
        denominatorX += (netWt * (edgeDegree - 1.0) *
                         _netsLinearizingWeightX[binIdx][binEdgeId]);
        denominatorY += (netWt * (edgeDegree - 1.0) *
                         _netsLinearizingWeightY[binIdx][binEdgeId]);
      } else {
        netXSum = netWt * netXSum;
        netYSum = netWt * netYSum;
        denominatorX += (netWt * (edgeDegree - 1.0));
        denominatorY += (netWt * (edgeDegree - 1.0));
      }
      xSum += netXSum;
      ySum += netYSum;
    }
  }

  if (denominatorX != 0 && denominatorY != 0) {
    double nodeSize =
        hgWDims.getNodeWidth(index) * hgWDims.getNodeHeight(index);
    xSum += nodeSize / _binNodesSizes[binIdx];
    ySum += nodeSize / _binNodesSizes[binIdx];
    xSum /= denominatorX;
    ySum /= denominatorY;
    optLoc.x = xSum;
    optLoc.y = ySum;
  } else
    optLoc = _nodesPlacementNatural[binIdx][_mapping[binIdx][index]];

  return optLoc;
}

/*TODO EXTEND THIS TO MULTI BIN CASE WITHOUT THAT IT WOULDN"T WORK*/
void AnalyticalSolver::computeLinearizingWeights() {
  const HGraphWDimensions& hgWDims = _rbplace.getHGraph();

  for (unsigned binIdx = 0; binIdx < _numBins; ++binIdx) {
    for (unsigned binNetId = 0; binNetId < _nets[binIdx].size(); ++binNetId) {
      unsigned netId = _nets[binIdx][binNetId];

      BBox netBBox;
      const HGFEdge& net = hgWDims.getEdgeByIdx(netId);

      itHGFNodeLocal nodeIt;
      unsigned nodeOffset;
      for (nodeIt = net.nodesBegin(), nodeOffset = 0; nodeIt != net.nodesEnd();
           nodeIt++, nodeOffset++) {
        const HGFNode& node = (**nodeIt);
        unsigned nodeIndex = node.getIndex();
        Point pOffset(0, 0);

        if (_usePinOffsets)
          pOffset = hgWDims.nodesOnEdgePinOffset(nodeOffset, netId,
                                                 _rbplace.getOrient(nodeIndex));

        if (_mapping[binIdx].count(nodeIndex) == 0)  // node not in this
                                                     // bin. so treat as fixed
        {
          double pinXLoc, pinYLoc;
          // node not in all bins , actual term
          if (_mappingNodesToBin.count(nodeIndex) == 0) {
            pinXLoc = _placement[nodeIndex].x + pOffset.x;
            pinYLoc = _placement[nodeIndex].y + pOffset.y;
          } else {
            unsigned nodeBinId = _mappingNodesToBin[nodeIndex];
            unsigned mappedIdx = _mapping[nodeBinId][nodeIndex];
            pinXLoc = _nodesPlacement[nodeBinId][mappedIdx].x + pOffset.x;
            pinYLoc = _nodesPlacement[nodeBinId][mappedIdx].y + pOffset.y;
          }
          // do terminal propagation
          if (pinXLoc < _binBBoxs[binIdx].xMin)
            pinXLoc = _binBBoxs[binIdx].xMin;
          else if (pinXLoc > _binBBoxs[binIdx].xMax)
            pinXLoc = _binBBoxs[binIdx].xMax;

          if (pinYLoc < _binBBoxs[binIdx].yMin)
            pinYLoc = _binBBoxs[binIdx].yMin;
          else if (pinYLoc > _binBBoxs[binIdx].yMax)
            pinYLoc = _binBBoxs[binIdx].yMax;

          netBBox.add(pinXLoc, pinYLoc);
        } else {
          unsigned mappedIdx = _mapping[binIdx][nodeIndex];
          if (_usePinOffsets) {
            pOffset = hgWDims.nodesOnEdgePinOffset(
                nodeOffset, netId, _rbplace.getOrient(nodeIndex));
          }

          netBBox.add(_nodesPlacement[binIdx][mappedIdx].x + pOffset.x,
                      _nodesPlacement[binIdx][mappedIdx].y + pOffset.y);
        }
      }
      if (!netBBox.isEmpty()) {
        _netsLinearizingWeightX[binIdx][binNetId] =
            std::max(netBBox.getWidth(), _avgCellSize);
        _netsLinearizingWeightY[binIdx][binNetId] =
            std::max(netBBox.getHeight(), _avgCellSize);

        _netsLinearizingWeightX[binIdx][binNetId] =
            1.0 / _netsLinearizingWeightX[binIdx][binNetId];
        _netsLinearizingWeightY[binIdx][binNetId] =
            1.0 / _netsLinearizingWeightY[binIdx][binNetId];
      }
    }
  }
}
