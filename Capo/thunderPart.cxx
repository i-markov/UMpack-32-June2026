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

#include "thunderPart.h"

#include <algorithm>
#include <deque>
#include <vector>

#include "Geoms/bbox.h"
#include "Geoms/point.h"
#include "PlaceEvals/pinPlEval.h"
#include "Stats/loadedDie.h"

using std::cout;
using std::deque;
using std::endl;
using std::string;
using std::vector;

const unsigned NO_BIN =
    UINT_MAX;  // special value to indicate node is not in a bin
const unsigned NOT_CUT = UINT_MAX;  // special value to indicate edge is not cut
const unsigned LOWEST_EDGE_DEGREE = 2;  // edge degrees
const unsigned NUM_EDGE_DEGREES = 4;    //  considered are [2,5]

ThunderPart::ThunderPart(vector<CapoBin*>& bins,
                         const HGraphWDimensions& hgraph,
                         const PlacementWOrient& placement,
                         const RBPlacement& rbplace, ThunderPartParams params)
    : _bins(bins),
      _rbplace(rbplace),
      _hgraph(hgraph),
      _placement(placement),
      _params(params) {
  initialize();

  if (_params.verb.getForActions() > 0)
    cout << "ThunderPart GO! (launching global placement refinement)" << endl;
}

void ThunderPart::partitionize(unsigned numMoves) {
  if (_params.verb.getForActions() > 0) {
    cout << "[";
    cout.flush();
  }

  Timer moveTimer;

  simpleGreedy(numMoves);
  // hillClimb(numMoves);

  moveTimer.stop();

  if (_params.verb.getForActions() > 0) {
    cout << "]" << endl;
    cout << numMoves << " ThunderMoves took: " << moveTimer << endl;
  }

  // update to put all the nodes in new bins, for the caller to use after
  // ThunderPart do two passes to preserve order of original nodes in bin
  vector<bool> seenNodes(_nodeBins.size(), false);
  for (unsigned b = 0; b < _bins.size(); b++)
    for (vector<unsigned>::const_iterator nodeIt = _bins[b]->cellIdsBegin();
         nodeIt != _bins[b]->cellIdsEnd(); nodeIt++)
      if (_nodeBins[*nodeIt] == b) {
        // still in same bin - push it to preserve order
        _resultBins[b].push_back(*nodeIt);
        seenNodes[*nodeIt] = true;
      }

  // push in all the nodes that have changed, in order seen
  for (unsigned nodeID = 0; nodeID < _nodeBins.size(); nodeID++)
    if (_nodeBins[nodeID] != NO_BIN && !seenNodes[nodeID])
      _resultBins[_nodeBins[nodeID]].push_back(nodeID);
}

void ThunderPart::recomputeExposedNets()
/*
 * look at all nets (edges), recompute them
 */
{
  for (itHGFEdgeGlobal edgeIt = _hgraph.edgesBegin();
       edgeIt != _hgraph.edgesEnd(); edgeIt++) {
    unsigned edgeIdx = (*edgeIt)->getIndex();
    updateExposedNet(edgeIdx);
  }
}

void ThunderPart::updateExposedNet(unsigned edgeIdx)
/*
 * look at one net (edgeIdx), update HPWL and netWeight.
 * Nets' weights are only 'exposed' if the nets are cut
 */
{
  updateHPWL(edgeIdx);

  const HGFEdge& edge = _hgraph.getEdgeByIdx(edgeIdx);
  unsigned edgeDegree = edge.getDegree();
  unsigned edgeSelectionBucketIdx;
  edgeSelectionBucketIdx = (edgeDegree < LOWEST_EDGE_DEGREE)
                               ? UINT_MAX
                               : edgeDegree - LOWEST_EDGE_DEGREE;

  // update data structures for selecting random edge
  if (netIsCut(edgeIdx)) {
    if (_edgeSelectionBucketLoc[edgeIdx] == NOT_CUT &&
        edgeSelectionBucketIdx >= 0 &&
        edgeSelectionBucketIdx < NUM_EDGE_DEGREES) {
      // edge wasn't cut before, and edge's degree is something we consider
      // - add it to the selection bucket
      _edgeSelectionBuckets[edgeSelectionBucketIdx].push_back(edgeIdx);
      _edgeSelectionBucketLoc[edgeIdx] =
          _edgeSelectionBuckets[edgeSelectionBucketIdx].size() - 1;
    }
  } else {
    if (!(_edgeSelectionBucketLoc[edgeIdx] == NOT_CUT)) {
      // this edge was cut before - remove it from the selection bucket
      unsigned bucketIdx = _edgeSelectionBucketLoc[edgeIdx];
      unsigned lastEdgeIdx =
          _edgeSelectionBuckets[edgeSelectionBucketIdx].size() - 1;

      // replace this edge with the last edge in the vector,
      // pop off the last element
      _edgeSelectionBuckets[edgeSelectionBucketIdx][bucketIdx] =
          _edgeSelectionBuckets[edgeSelectionBucketIdx][lastEdgeIdx];
      _edgeSelectionBuckets[edgeSelectionBucketIdx].pop_back();
      _edgeSelectionBucketLoc[edgeIdx] = NOT_CUT;
    }
  }
}

double ThunderPart::getHPWL() { return _totalHPWL; }

void ThunderPart::updateHPWL(unsigned edgeIdx)
/*
 * incrementally update HPWL by re-evaluating
 * the bbox of nodes on this edge
 */
{
  _totalHPWL -= _HPWL[edgeIdx];

  bool usePinLocs = true;
  bool useWts = true;
  _HPWL[edgeIdx] =
      pinPlEval::oneNetHPWL(_placement, _hgraph, edgeIdx, usePinLocs, useWts);

  _totalHPWL += _HPWL[edgeIdx];
}

bool ThunderPart::netIsCut(unsigned edgeIdx) const
/*
 * if 2 nodes on an edge are in different bins, return true
 */
{
  const HGFEdge& edge = _hgraph.getEdgeByIdx(edgeIdx);

  unsigned binID = UINT_MAX;

  itHGFNodeLocal node;
  for (node = edge.nodesBegin(); node != edge.nodesEnd(); node++) {
    unsigned nodeID = (*node)->getIndex();

    if (_nodeBins[nodeID] == NO_BIN)
      continue;
    else if (binID == UINT_MAX) {  // initialization
      binID = _nodeBins[nodeID];
      continue;
    } else if (binID != _nodeBins[nodeID])
      return true;
  }

  return false;
}

void ThunderPart::initialize() {
  unsigned numNets = _hgraph.getNumEdges();

  _HPWL.resize(numNets);
  fill(_HPWL.begin(), _HPWL.end(), 0.0);
  _resultBins.resize(_bins.size());

  _edgeSelectionBuckets.resize(NUM_EDGE_DEGREES);
  _edgeSelectionBucketLoc.resize(numNets);
  fill(_edgeSelectionBucketLoc.begin(), _edgeSelectionBucketLoc.end(), NOT_CUT);

  _totalHPWL = 0.0;

  _nodeBins.resize(_hgraph.getNumNodes());
  fill(_nodeBins.begin(), _nodeBins.end(), NO_BIN);
  _binCellAreas.resize(_bins.size());
  fill(_binCellAreas.begin(), _binCellAreas.end(), 0.0);

  for (unsigned b = 0; b < _bins.size(); b++)
    if (_bins[b] != NULL) {
      for (vector<unsigned>::const_iterator nodeIt = _bins[b]->cellIdsBegin();
           nodeIt != _bins[b]->cellIdsEnd(); nodeIt++) {
        _nodeBins[*nodeIt] = b;
      }

      _binCellAreas[b] = _bins[b]->getTotalCellArea();
    }

  recomputeExposedNets();
}

void ThunderPart::moveCell(unsigned nodeID, unsigned binID, Point newLoc)
/*
 * move nodeID to binID
 */
{
  if (isOutOfCore(newLoc)) {
    cout << "FATAL: cell " << nodeID << " wants to go to " << newLoc << endl;
    cout << "dest bin bbox: " << _bins[binID]->getBBox();
    cout << "core bbox " << _rbplace.getBBox(false);
    abkwarn(0, "cell is being moved out of core");
  }

  // update bin assignment
  unsigned oldBinID = _nodeBins[nodeID];
  _nodeBins[nodeID] = binID;

  // reassign location
  _placement[nodeID] = newLoc;

  // update book-keeping

  double cellArea =
      _hgraph.getNodeWidth(nodeID) * _hgraph.getNodeHeight(nodeID);
  _binCellAreas[oldBinID] -= cellArea;
  _binCellAreas[binID] += cellArea;

  HGFNode node = _hgraph.getNodeByIdx(nodeID);
  for (itHGFEdgeLocal adjEdgesIt = node.edgesBegin();
       adjEdgesIt != node.edgesEnd(); ++adjEdgesIt) {
    updateExposedNet((*adjEdgesIt)->getIndex());
  }
}

void ThunderPart::simpleGreedy(unsigned numMoves) {
  double oldHPWL = getHPWL();

  for (unsigned moveCount = 0; moveCount < numMoves; ++moveCount) {
    if (noVisibleEdges()) continue;

    const HGFEdge& edge = selectVisibleEdge();

    deque<unsigned> movableNodesOnEdge;
    vector<unsigned> binsSpannedByEdge;

    // find candidate nodes and bins for moves
    BBox edgeBBox;
    for (itHGFNodeLocal nodeIt = edge.nodesBegin(); nodeIt != edge.nodesEnd();
         nodeIt++) {
      unsigned nodeID = (*nodeIt)->getIndex();
      if (_nodeBins[nodeID] == NO_BIN) {
        edgeBBox += _placement[nodeID];
        continue;
      }

      if (isMovable(nodeID)) {
        if (!edgeBBox.hasInside(_placement[nodeID])) {
          // this node expands the bbox of this edge
          // - give it higher priority
          movableNodesOnEdge.push_front(nodeID);
        } else {
          // don't really care about this node's priority
          movableNodesOnEdge.push_back(nodeID);
        }
      }

      if (!_hgraph.isTerminal(nodeID))
        binsSpannedByEdge.push_back(_nodeBins[nodeID]);

      edgeBBox += _placement[nodeID];

      // find adjacent bins
      const HGFNode& node = _hgraph.getNodeByIdx(nodeID);
      for (itHGFEdgeLocal edgeIt = node.edgesBegin(); edgeIt != node.edgesEnd();
           edgeIt++)
        for (itHGFNodeLocal nodeIt = (*edgeIt)->nodesBegin();
             nodeIt != (*edgeIt)->nodesEnd(); nodeIt++) {
          unsigned neighborIdx = (*nodeIt)->getIndex();
          if (_nodeBins[neighborIdx] == NO_BIN) continue;

          // push adjacent bin. NOTE: repeated bins are intentionally
          // allowed. The intuition is that repeated bins are more
          // likely to be selected, and this is desired because repeated
          // bins are bins which have multiple nodes, and so it seems
          // that we'd want to move distant nodes to where most of the
          // other nodes are.
          if (!_hgraph.isTerminal(neighborIdx))
            binsSpannedByEdge.push_back(_nodeBins[neighborIdx]);
        }
    }

    if (movableNodesOnEdge.size() < 2) continue;
    if (binsSpannedByEdge.size() < 2) continue;

    // probabilistically pick a node for moving, prefer nodes at front
    // (these are the ones that expanded the edge's bbox during construction)
    unsigned selectedNode;
    for (;;)
      if (coinFlip()) {
        selectedNode = movableNodesOnEdge.front();
        break;
      } else {
        unsigned unpicked = movableNodesOnEdge.front();
        movableNodesOnEdge.pop_front();
        movableNodesOnEdge.push_back(unpicked);
      }

    unsigned binMovedFrom = _nodeBins[selectedNode];
    Point pointMovedFrom = _placement[selectedNode];

    // try find the bestest bin for this node using weber,
    // but only use this sometimes
    unsigned weberBin;
    if (coinFlip() && (weberBin = getWeberBin(selectedNode)) != NO_BIN &&
        weberBin != _nodeBins[selectedNode] &&
        !willOverfillBin(selectedNode, weberBin)) {
      // found the most awesomest bin for this node,
      // wrt this edge

      double oldAvgBinTargetUtilDiff =
          ((_binCellAreas[binMovedFrom] -
            _bins[binMovedFrom]->getTotalCellArea()) +
           (_binCellAreas[weberBin] - _bins[weberBin]->getTotalCellArea())) /
          2.0;

      moveCell(selectedNode, weberBin,
               findBestCaseDest(selectedNode, weberBin));

      double newHPWL = getHPWL();

      double newAvgBinTargetUtilDiff =
          ((_binCellAreas[binMovedFrom] -
            _bins[binMovedFrom]->getTotalCellArea()) +
           (_binCellAreas[weberBin] - _bins[weberBin]->getTotalCellArea())) /
          2.0;

      if (newHPWL < oldHPWL ||
          (newHPWL == oldHPWL &&
           newAvgBinTargetUtilDiff < oldAvgBinTargetUtilDiff)) {
        if (newHPWL < oldHPWL) announceHPWL('*', newHPWL);

        oldHPWL = newHPWL;

        continue;  // next move!
      } else {
        // move may be good for this edge, but bad globally - roll back
        moveCell(selectedNode, binMovedFrom, pointMovedFrom);
      }
      /* FALLTHROUGH */
    }

    unsigned randPick =
        unsigned(binsSpannedByEdge.size() / (UINT_MAX + 1.0) * _randomNum);
    unsigned selectedBin = binsSpannedByEdge[randPick];

    unsigned binMovedTo = selectedBin;
    Point pointMovedTo = findBestCaseDest(selectedNode, selectedBin);

    if (willOverfillBin(selectedNode, selectedBin)) continue;

    double oldAvgBinTargetUtilDiff =
        ((_binCellAreas[binMovedFrom] -
          _bins[binMovedFrom]->getTotalCellArea()) +
         (_binCellAreas[binMovedTo] - _bins[binMovedTo]->getTotalCellArea())) /
        2.0;

    moveCell(selectedNode, selectedBin, pointMovedTo);

    double newHPWL = getHPWL();

    double newAvgBinTargetUtilDiff =
        ((_binCellAreas[binMovedFrom] -
          _bins[binMovedFrom]->getTotalCellArea()) +
         (_binCellAreas[binMovedTo] - _bins[binMovedTo]->getTotalCellArea())) /
        2.0;

    // keep move if new wirelength is better, or if wirelength stays the
    // same but is closer to Capo's original bin utilization
    if (newHPWL < oldHPWL ||
        (newHPWL == oldHPWL &&
         newAvgBinTargetUtilDiff < oldAvgBinTargetUtilDiff)) {
      // good move, keep
      if (newHPWL < oldHPWL) announceHPWL('.', newHPWL);

      oldHPWL = newHPWL;
    } else {
      // bad move - roll back
      moveCell(selectedNode, binMovedFrom, pointMovedFrom);
    }
  }
}

#if 0
void
ThunderPart::hillClimb(unsigned numMoves) 
{
	PlacementWOrient bestPlacement(_placement);
	vector<int> bestNodeBins(_nodeBins);
	double minHPWL = getHPWL();

	typedef enum { UP, DOWN } Direction; // hill-climbing direction
	Direction direction = DOWN; // start downward initially

	double oldHPWL = minHPWL;
	
	for (unsigned moveCount = 0; moveCount < numMoves; ++moveCount) {
		unsigned lowestBinMovedTo;
		Point lowestPointMovedTo;
		unsigned lowestNodeIdx;
		double lowestHPWL = DBL_MAX;
		bool lowMoveSet = false;

		unsigned highestBinMovedTo;
		Point highestPointMovedTo;
		unsigned highestNodeIdx;
		double highestHPWL = -DBL_MAX;
		bool highMoveSet = false;

		// check a few random edges
		for (unsigned numEdgestoCheck = 0; numEdgestoCheck < 5; numEdgestoCheck++) {
			if (noVisibleEdges())
				continue;

			const HGFEdge& edge = selectVisibleEdge();
			
			// pick a random node on the edge to move

			vector<unsigned> nodesOnEdge;
			for (itHGFNodeLocal nodeIt = edge.nodesBegin();
				nodeIt != edge.nodesEnd(); nodeIt++)
				nodesOnEdge.push_back((*nodeIt)->getIndex());

			random_shuffle(nodesOnEdge.begin(), nodesOnEdge.end(), _randomNum);

			unsigned randNode;
			unsigned nodeIdx;
			for (randNode = 0; randNode < nodesOnEdge.size(); randNode++) {
				unsigned nodeID = nodesOnEdge[randNode];
				if (_nodeBins[nodeID] == NO_BIN)
					continue;
				if (!_hgraph.isNodeMacro(nodeID)) {
					continue;

				nodeIdx = nodeID;
				break;
			}
			if (randNode == nodesOnEdge.size()) {
				// no valid cells found for moving
				//abkwarn(0, "aaronnn: no node found to act upon");
				continue;
			}

			unsigned binMovedFrom = 0;
			unsigned binMovedTo = 0;
			Point pointMovedFrom;
			Point pointMovedTo;

			// consider all the bins on the selected edge for hill-climbing
			itHGFNodeLocal nodeIt;
			for (nodeIt = edge.nodesBegin();
				nodeIt != edge.nodesEnd(); ++nodeIt)
			{
				if (_nodeBins[(*nodeIt)->getIndex()] == NO_BIN)
					continue;
				else if (_nodeBins[(*nodeIt)->getIndex()] == _nodeBins[nodeIdx])
					continue;
				else
				{
					// find hill-climbing moves
					binMovedFrom = _nodeBins[nodeIdx];
					binMovedTo = _nodeBins[(*nodeIt)->getIndex()];
					pointMovedFrom = _placement[nodeIdx];
					pointMovedTo = findBestCaseDest(nodeIdx, binMovedTo);

					moveCell(nodeIdx, binMovedFrom, pointMovedTo);

					double newHPWL = getHPWL();

					double capacityUsed = _binCellAreas[binMovedTo];
					if (capacityUsed > _bins[binMovedTo]->getMaxRecommendedCellArea())
					{
						// area constraints violated
						// don't save this move
					} else {
						// remember the lowest highest-cost move
						if ((!highMoveSet && newHPWL > oldHPWL) 
							|| (newHPWL > oldHPWL && newHPWL < highestHPWL)) {
							highestBinMovedTo = binMovedTo;
							highestPointMovedTo = pointMovedTo;
							highestNodeIdx = nodeIdx;
							highestHPWL = newHPWL;
							highMoveSet = true;
						} 

						// remember least-cost move
						if (newHPWL < lowestHPWL) {
							lowestBinMovedTo = binMovedTo;
							lowestPointMovedTo = pointMovedTo;
							lowestNodeIdx = nodeIdx;
							lowestHPWL = newHPWL;
							lowMoveSet = true;
						}
					}

					// move cell back
					moveCell(nodeIdx, binMovedFrom, pointMovedFrom);
				}
			}
		}

		if (!lowMoveSet || !highMoveSet) {
			// not enough hill-climbing moves found
			//abkwarn(0, "aaronnn: no hill-climbing moves found");
			continue;
		}

		double newHPWL = oldHPWL;

		// decide which way to climb this hill, and DOo eeeeeT!
		if (direction == DOWN) {
			if (lowestHPWL <= oldHPWL) {
				// keep going down
				moveCell(lowestNodeIdx, lowestBinMovedTo, lowestPointMovedTo);
				newHPWL = lowestHPWL;
			} else {
				// start going up
				moveCell(highestNodeIdx, highestBinMovedTo, highestPointMovedTo);
				newHPWL = highestHPWL;
				//direction = UP;
				//cout << "DBG going up at " << oldHPWL << endl;
			}
		} else { // direction == UP
			if (highestHPWL > oldHPWL) {
				// keep going up
				moveCell(highestNodeIdx, highestBinMovedTo, highestPointMovedTo);
				newHPWL = highestHPWL;
			} else {
				// start going down
				moveCell(lowestNodeIdx, lowestBinMovedTo, lowestPointMovedTo);
				newHPWL = lowestHPWL;
				direction = DOWN;
				//cout << "DBG going down at " << oldHPWL << endl;
			}
		}

		if (newHPWL < minHPWL) {
			announceHPWL('.', newHPWL);

			// remember best placement
			minHPWL = newHPWL;
			bestPlacement = _placement;
			bestNodeBins = _nodeBins;
			// this makes a copy of the placement.
			// TODO remember moves made instead?
		}

		oldHPWL = newHPWL;
	}

	// update the active placement
	_placement = bestPlacement;
	_nodeBins = bestNodeBins;
	recomputeExposedNets();
}
#endif

const HGFEdge& ThunderPart::selectVisibleEdge() {
  // first, find the edge degree we want

  // get bucket weights
  std::vector<double> edgeBucketWeights(NUM_EDGE_DEGREES);
  for (unsigned i = 0; i < NUM_EDGE_DEGREES; i++) {
    double netWeight;

    // manually adjust weights
    if (i ==
        0)  // ... because Capo's bisection would have already done a good job:
      netWeight = 1000 / (1 + LOWEST_EDGE_DEGREE - 1);
    else if (i == 1)  // ... because they potentially offer easiest improvement:
      netWeight = 1000 / (0 + LOWEST_EDGE_DEGREE - 1);
    else if (i == 2)  // ... because they are harder to improve:
      netWeight = 1000 / (2 + LOWEST_EDGE_DEGREE - 1);
    else if (i == 3)  // ... because they are harder to improve:
      netWeight = 1000 / (3 + LOWEST_EDGE_DEGREE - 1);
    else
      netWeight = 1000 / (i + LOWEST_EDGE_DEGREE - 1);

    edgeBucketWeights[i] = _edgeSelectionBuckets[i].size() * netWeight;
  }

  // pick edge bucket
  unsigned selectedEdgeBucket = LoadedDie(edgeBucketWeights, _randomNum);

  // finally, find a cut edge that has this edge degree

  unsigned selectedEdge =
      unsigned(_edgeSelectionBuckets[selectedEdgeBucket].size() /
               (UINT_MAX + 1.0) * _randomNum);

  return _hgraph.getEdgeByIdx(
      _edgeSelectionBuckets[selectedEdgeBucket][selectedEdge]);
}

Point ThunderPart::findBestCaseDest(unsigned nodeIdx, unsigned binMovedTo) {
  return _bins[binMovedTo]->getCenter();
  /*
          unsigned binMovedFrom = _nodeBins[nodeIdx];
          Point pointMovedFrom = _placement[nodeIdx];
          double bestHPWL = DBL_MAX;
          Point pointMovedTo;
          Point bestCasePointMovedTo;
          double HPWL;

          pointMovedTo = _bins[binMovedTo]->getCenter();
          moveCell(nodeIdx, binMovedTo, pointMovedTo);
          HPWL = getHPWL();
          if (HPWL < bestHPWL) {
                  bestCasePointMovedTo = pointMovedTo;
                  bestHPWL = HPWL;
          }
          moveCell(nodeIdx, binMovedFrom, pointMovedFrom);

          BBox binBBox = _bins[binMovedTo]->getBBox();
          pointMovedTo = Point(binBBox.xMin, binBBox.yMin);
          moveCell(nodeIdx, binMovedTo, pointMovedTo);
          HPWL = getHPWL();
          if (HPWL < bestHPWL) {
                  bestCasePointMovedTo = pointMovedTo;
                  bestHPWL = HPWL;
          }
          moveCell(nodeIdx, binMovedFrom, pointMovedFrom);

          pointMovedTo = Point(binBBox.xMax, binBBox.yMax);
          moveCell(nodeIdx, binMovedTo, pointMovedTo);
          HPWL = getHPWL();
          if (HPWL < bestHPWL) {
                  bestCasePointMovedTo = pointMovedTo;
                  bestHPWL = HPWL;
          }
          moveCell(nodeIdx, binMovedFrom, pointMovedFrom);

          pointMovedTo = Point(binBBox.xMin, binBBox.yMax);
          moveCell(nodeIdx, binMovedTo, pointMovedTo);
          HPWL = getHPWL();
          if (HPWL < bestHPWL) {
                  bestCasePointMovedTo = pointMovedTo;
                  bestHPWL = HPWL;
          }
          moveCell(nodeIdx, binMovedFrom, pointMovedFrom);

          pointMovedTo = Point(binBBox.xMax, binBBox.yMin);
          moveCell(nodeIdx, binMovedTo, pointMovedTo);
          HPWL = getHPWL();
          if (HPWL < bestHPWL) {
                  bestCasePointMovedTo = pointMovedTo;
                  bestHPWL = HPWL;
          }
          moveCell(nodeIdx, binMovedFrom, pointMovedFrom);

          return bestCasePointMovedTo;
  */
}

bool ThunderPart::noVisibleEdges() {
  for (unsigned i = 0; i < NUM_EDGE_DEGREES; i++)
    if (_edgeSelectionBuckets[i].size())  // we have cut/visible edges
      return false;

  return true;
}

unsigned ThunderPart::getWeberBin(unsigned nodeIdx)
/*
 * look at all edges connected to this node,
 * relocate node to where it wants to go
 * returns the bin that nodeIdx is currently in
 * if no destination bin is found.
 */
{
  vector<bool> visited(_hgraph.getNumEdges(), false);  // incident to bin
  vector<unsigned> adjacentBins;

  const HGFNode& node = _hgraph.getNodeByIdx(nodeIdx);

  vector<double> xEnds, yEnds;
  xEnds.reserve(2 * node.getDegree()), yEnds.reserve(2 * node.getDegree());
  for (itHGFEdgeLocal edgeIt = node.edgesBegin(); edgeIt != node.edgesEnd();
       edgeIt++) {
    BBox box;  // created empty
    const HGFEdge& edge = (**edgeIt);
    unsigned edgeIdx = edge.getIndex();

    if (visited[edgeIdx]) continue;

    visited[edgeIdx] = true;
    bool hasFixed = false;

    for (itHGFNodeLocal nodeIt = edge.nodesBegin(); nodeIt != edge.nodesEnd();
         nodeIt++) {
      unsigned nodeNum = (*nodeIt)->getIndex();

      if (_nodeBins[nodeNum] == NO_BIN) {
        // we can't move the cell to this cell's bin,
        // since this cell isn't in a bin
        continue;
      }

      if (nodeNum != nodeIdx) {
        hasFixed = true;
        box += _placement[nodeNum];
        adjacentBins.push_back(_nodeBins[nodeNum]);
      }
    }

    if (hasFixed) {
      xEnds.push_back(box.xMin);
      xEnds.push_back(box.xMax);
      yEnds.push_back(box.yMin);
      yEnds.push_back(box.yMax);
    }
  }

  if (xEnds.empty()) return _nodeBins[nodeIdx];

  sort(xEnds.begin(), xEnds.end());
  sort(yEnds.begin(), yEnds.end());

  unsigned n = xEnds.size() / 2;
  double xLeft = xEnds[n - 1], xRight = xEnds[n];
  double yLeft = yEnds[n - 1], yRight = yEnds[n];

  BBox weberBin(xLeft, yLeft, xRight, yRight);

  for (unsigned i = 0; i < adjacentBins.size(); i++) {
    unsigned binIdx = adjacentBins[i];
    if (_bins[binIdx]->getBBox().intersects(weberBin)) {
      return binIdx;
    }
  }

  return _nodeBins[nodeIdx];
}

bool ThunderPart::willOverfillBin(unsigned nodeID, unsigned binID) {
  if (binID == NO_BIN) {
    // you can never overfill an imaginary bin
    return false;
  }

  double cellArea =
      _hgraph.getNodeWidth(nodeID) * _hgraph.getNodeHeight(nodeID);
  double expectedCapacityUsed = _binCellAreas[binID] + cellArea;

  if (expectedCapacityUsed > _bins[binID]->getMaxRecommendedCellArea())
    return true;

  return false;
}

bool ThunderPart::coinFlip() {
  if (double(_randomNum) / double(UINT_MAX) > 0.5) return true;

  return false;
}

bool ThunderPart::coinFlip(double sideProb) {
  if (double(_randomNum) / double(UINT_MAX) > sideProb) return true;

  return false;
}

void ThunderPart::announceHPWL(char c, double HPWL) {
  if (_params.verb.getForActions() > 0) {
    cout << c;
    cout.flush();
    // cout << "########!!! T-H-U-N-D-E-R !!!######## (" << HPWL << ")" << endl;
  }
}

bool ThunderPart::isOutOfCore(Point& point) {
  bool includeTerminals = false;
  if (!_rbplace.getBBox(includeTerminals).contains(point)) return true;

  return false;
}

bool ThunderPart::isMovable(unsigned nodeID) {
  if (!_hgraph.isNodeMacro(nodeID) && !_hgraph.isTerminal(nodeID) &&
      !_rbplace.isFixed(nodeID))
    return true;

  return false;
}
