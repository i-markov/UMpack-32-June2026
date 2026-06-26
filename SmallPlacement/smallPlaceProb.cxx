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

//  Created by Igor Markov, April 30, 1999
//  CHANGES

#include "smallPlaceProb.h"

#include <unordered_map>
#include "HGraphWDims/hgWDims.h"
#include "Placement/placeWOri.h"
#include "baseSmallPlPr.h"

using std::vector;

SmallPlacementProblem::SmallPlacementProblem(
    const HGraphWDimensions& hgwd, const PlacementWOrient& pl,
    const bit_vector& placed, const vector<unsigned>& movables,
    const SmallPlacementRow& row, SmallPlacementBitBoardContainer& bitBoards,
    Verbosity verb)

    : BaseSmallPlacementProblem(movables.size(), verb), _bitBoards(bitBoards) {
  _fixed = bit_vector(movables.size(), 0);
  _rows.push_back(row);

  unsigned k, numCells = movables.size();
  _bitBoards.mobileCells.reset(hgwd.getNumNodes());
  for (k = 0; k != numCells; k++) _bitBoards.mobileCells.setBit(movables[k]);

  _cellWidths.reserve(numCells);
  for (k = 0; k != numCells; k++) {
    _cellWidths.push_back(hgwd.getNodeWidth(movables[k]));
    abkassert(_cellWidths[k] != -DBL_MAX, "uninitialized cell width");
  }

  _bitBoards.seenNets.reset(hgwd.getNumEdges());
  _bitBoards.essentialNets.reset(hgwd.getNumEdges());

  for (k = 0; k != numCells; k++) {
    const HGFNode& node = hgwd.getNodeByIdx(movables[k]);
    for (itHGFEdgeLocal e = node.edgesBegin(); e != node.edgesEnd(); e++) {
      const HGFEdge& edge = (**e);
      unsigned netIdx = edge.getIndex();
      if (!_bitBoards.seenNets.isBitSet(netIdx)) {
        _bitBoards.seenNets.setBit(netIdx);
        unsigned haveFixed = false;
        BBox box;
        unsigned nodeOffset = 0;
        for (itHGFNodeLocal otherNode = edge.nodesBegin();
             otherNode != edge.nodesEnd(); otherNode++, nodeOffset++) {
          unsigned nodeIndex = (*otherNode)->getIndex();
          if (_bitBoards.mobileCells.isBitSet(nodeIndex)) continue;

          if (placed[nodeIndex]) {
            //                    Point pinOffset=hgwd.nodesOnEdgePinOffset
            //				(nodeOffset,netIdx,pl.getOrient(nodeIndex));
            //                    box+=(pl[nodeIndex]+pinOffset);
            BBox pinOffset = hgwd.nodesOnEdgePinOffsetBBox(
                nodeOffset, netIdx, pl.getOrient(nodeIndex));
            box += (pl[nodeIndex] + pinOffset.getBottomLeft());
            box += (pl[nodeIndex] + pinOffset.getTopRight());
          } else
            box += pl[nodeIndex];

          haveFixed = true;
        }

        if (!haveFixed || box.xMin > row.xMin || box.xMax < row.xMax) {
          _netTerminalBBoxes.push_back(box);
          _bitBoards.essentialNets.setBit(netIdx);
        }
      }
    }
  }

  const auto& essNetIds = _bitBoards.essentialNets.getIndicesOfSetBits();
  unsigned numNets = essNetIds.size();

  _netlist.setup(numCells, numNets);

  std::unordered_map<unsigned, unsigned, std::hash<unsigned>, std::equal_to<unsigned> >
      origToNewIdx;
  for (k = 0; k != numCells; k++) origToNewIdx[movables[k]] = k;

  // setup the pin offsets for the movable cells
  for (k = 0; k != numNets; k++) {
    unsigned netIndex = essNetIds[k];
    const HGFEdge& edge = hgwd.getEdgeByIdx(netIndex);
    unsigned nodeOffset = 0;
    for (itHGFNodeLocal node = edge.nodesBegin(); node != edge.nodesEnd();
         node++, nodeOffset++) {
      unsigned nodeIndex = (*node)->getIndex();
      if (!_bitBoards.mobileCells.isBitSet(nodeIndex)) continue;
      //          Point pinOffset=hgwd.nodesOnEdgePinOffset
      //				(nodeOffset,netIndex,pl.getOrient(nodeIndex));
      BBox pinOffset = hgwd.nodesOnEdgePinOffsetBBox(nodeOffset, netIndex,
                                                     pl.getOrient(nodeIndex));

      _netlist.setPinOffset(origToNewIdx[nodeIndex], k, pinOffset);
    }
  }

  _bitBoards.seenNets.clear();
  _bitBoards.essentialNets.clear();
  _bitBoards.mobileCells.clear();
}

// [MAR] get cell width snapped to fill site spacing (not site width)
// so that cells will always be placed on site boundaries.
double SmallPlacementProblem::getSnappedCellWidth(const HGraphWDimensions& hgwd,
                                                  unsigned k) const {
  double width = hgwd.getNodeWidth(k);
  double spacing = _rows[0].getSiteInterval();
  double nSites = ceil(width / spacing);
  return nSites * spacing;
}

// this constructor is added by saurabh adya to take care of fake cells(WS)
// in the smallplacement problem
//
//  [MAR] this is the constructor used in Capo refineBin.cxx
//

SmallPlacementProblem::SmallPlacementProblem(
    const HGraphWDimensions& hgwd, const PlacementWOrient& pl,
    const bit_vector& placed, const bit_vector& fixed,
    const vector<unsigned>& movables, const vector<SmallPlacementRow>& rows,
    const vector<double>& whiteSpaceWidths,
    SmallPlacementBitBoardContainer& bitBoards, Verbosity verb)

    : BaseSmallPlacementProblem(movables.size(), verb), _bitBoards(bitBoards) {
  _fixed = fixed;
  unsigned i;
  for (i = 0; i < rows.size(); ++i) _rows.push_back(rows[i]);

  unsigned k, numCells = movables.size();
  _bitBoards.mobileCells.reset(hgwd.getNumNodes());

  for (k = 0; k != numCells; k++)
    if (movables[k] != UINT_MAX) _bitBoards.mobileCells.setBit(movables[k]);

  double rowsXMax = -DBL_MAX;
  double rowsXMin = DBL_MAX;
  double rowsYMax = -DBL_MAX;
  double rowsYMin = DBL_MAX;

  for (i = 0; i < rows.size(); ++i) {
    if (rows[i].xMax > rowsXMax) rowsXMax = rows[i].xMax;
    if (rows[i].xMin < rowsXMin) rowsXMin = rows[i].xMin;

    if (rows[i].yMax > rowsYMax) rowsYMax = rows[i].yMax;
    if (rows[i].yMin < rowsYMin) rowsYMin = rows[i].yMin;
  }

  _cellWidths.reserve(numCells);
  unsigned whiteSpaceCtr = 0;
  for (k = 0; k != numCells; k++) {
    if (movables[k] != UINT_MAX) {
      _cellWidths.push_back(getSnappedCellWidth(hgwd, movables[k]));
    } else {
      _cellWidths.push_back(whiteSpaceWidths[whiteSpaceCtr]);
      ++whiteSpaceCtr;
    }
    abkassert(_cellWidths[k] != -DBL_MAX, "uninitialized cell width");
  }

  _bitBoards.seenNets.reset(hgwd.getNumEdges());
  _bitBoards.essentialNets.reset(hgwd.getNumEdges());

  for (k = 0; k != numCells; k++) {
    if (movables[k] != UINT_MAX && _fixed[k] != 1) {
      const HGFNode& node = hgwd.getNodeByIdx(movables[k]);
      for (itHGFEdgeLocal e = node.edgesBegin(); e != node.edgesEnd(); e++) {
        const HGFEdge& edge = (**e);
        unsigned netIdx = edge.getIndex();
        if (!_bitBoards.seenNets.isBitSet(netIdx)) {
          _bitBoards.seenNets.setBit(netIdx);
          unsigned haveFixed = false;
          BBox box;
          unsigned nodeOffset = 0;
          for (itHGFNodeLocal otherNode = edge.nodesBegin();
               otherNode != edge.nodesEnd(); otherNode++, nodeOffset++) {
            unsigned nodeIndex = (*otherNode)->getIndex();
            if (_bitBoards.mobileCells.isBitSet(nodeIndex)) continue;

            if (placed[nodeIndex]) {
              //			 Point
              //pinOffset=hgwd.nodesOnEdgePinOffset
              //			   (nodeOffset,netIdx,pl.getOrient(nodeIndex));
              //			 box+=(pl[nodeIndex]+pinOffset);
              BBox pinOffset = hgwd.nodesOnEdgePinOffsetBBox(
                  nodeOffset, netIdx, pl.getOrient(nodeIndex));
              box += (pl[nodeIndex] + pinOffset.getBottomLeft());
              box += (pl[nodeIndex] + pinOffset.getTopRight());
            } else {
              box += pl[nodeIndex];
            }
            haveFixed = true;
          }
          // remove inessential nets
          if (!haveFixed || (box.xMin > rowsXMin || box.xMax < rowsXMax ||
                             box.yMin > rowsYMin || box.yMax < rowsYMax))

          {
            _netTerminalBBoxes.push_back(box);
            _bitBoards.essentialNets.setBit(netIdx);
          }
        }
      }
    }
  }
  const auto& essNetIds = _bitBoards.essentialNets.getIndicesOfSetBits();
  unsigned numNets = essNetIds.size();

  _netlist.setup(numCells, numNets);

  std::unordered_map<unsigned, unsigned, std::hash<unsigned>, std::equal_to<unsigned> >
      origToNewIdx;
  for (k = 0; k != numCells; k++)
    if (movables[k] != UINT_MAX) origToNewIdx[movables[k]] = k;

  // setup the pin offsets for the movable cells
  for (k = 0; k != numNets; k++) {
    unsigned netIndex = essNetIds[k];
    const HGFEdge& edge = hgwd.getEdgeByIdx(netIndex);
    unsigned nodeOffset = 0;
    for (itHGFNodeLocal node = edge.nodesBegin(); node != edge.nodesEnd();
         node++, nodeOffset++) {
      unsigned nodeIndex = (*node)->getIndex();
      if (!_bitBoards.mobileCells.isBitSet(nodeIndex)) continue;
      //          Point pinOffset=hgwd.nodesOnEdgePinOffset
      //				(nodeOffset,netIndex,pl.getOrient(nodeIndex));

      BBox pinOffset = hgwd.nodesOnEdgePinOffsetBBox(nodeOffset, netIndex,
                                                     pl.getOrient(nodeIndex));

      _netlist.setPinOffset(origToNewIdx[nodeIndex], k, pinOffset);
    }
  }

  _bitBoards.seenNets.clear();
  _bitBoards.essentialNets.clear();
  _bitBoards.mobileCells.clear();
}

// this constructor is added by saurabh adya to allow clustering of cells
// in the smallplacement problem

SmallPlacementProblem::SmallPlacementProblem(
    const HGraphWDimensions& hgwd, const PlacementWOrient& pl,
    const bit_vector& placed, const bit_vector& fixed,
    const vector<vector<unsigned> >& movables,
    const vector<SmallPlacementRow>& rows,
    const vector<double>& whiteSpaceWidths,
    SmallPlacementBitBoardContainer& bitBoards, Verbosity verb)

    : BaseSmallPlacementProblem(movables.size(), verb), _bitBoards(bitBoards) {
  _fixed = fixed;
  unsigned i;
  for (i = 0; i < rows.size(); ++i) _rows.push_back(rows[i]);

  unsigned k, numCells = movables.size();
  _bitBoards.mobileCells.reset(hgwd.getNumNodes());
  _bitBoards.clusteredCells.reset(hgwd.getNumNodes());

  for (k = 0; k != numCells; k++)
    if (movables[k][0] != UINT_MAX)  // if not whitespace
    {
      for (i = 0; i < movables[k].size(); ++i) {
        if (movables[k][i] != UINT_MAX)
          _bitBoards.mobileCells.setBit(movables[k][i]);
        if (movables[k].size() > 1)
          _bitBoards.clusteredCells.setBit(movables[k][i]);
      }
    }

  double rowsXMax = -DBL_MAX;
  double rowsXMin = DBL_MAX;
  double rowsYMax = -DBL_MAX;
  double rowsYMin = DBL_MAX;

  for (i = 0; i < rows.size(); ++i) {
    if (rows[i].xMax > rowsXMax) rowsXMax = rows[i].xMax;
    if (rows[i].xMin < rowsXMin) rowsXMin = rows[i].xMin;

    if (rows[i].yMax > rowsYMax) rowsYMax = rows[i].yMax;
    if (rows[i].yMin < rowsYMin) rowsYMin = rows[i].yMin;
  }

  _cellWidths.reserve(numCells);
  unsigned whiteSpaceCtr = 0;
  double nodeWidth = 0;

  for (k = 0; k != numCells; k++) {
    if (movables[k][0] != UINT_MAX)  // if not whitespace
    {
      nodeWidth = 0;
      for (i = 0; i < movables[k].size(); ++i)
        nodeWidth += hgwd.getNodeWidth(movables[k][i]);

      _cellWidths.push_back(nodeWidth);
    } else {
      _cellWidths.push_back(whiteSpaceWidths[whiteSpaceCtr]);
      ++whiteSpaceCtr;
    }
    abkassert(_cellWidths[k] != -DBL_MAX, "uninitialized cell width");
  }

  _bitBoards.seenNets.reset(hgwd.getNumEdges());
  _bitBoards.essentialNets.reset(hgwd.getNumEdges());

  for (k = 0; k != numCells; k++) {
    if (movables[k][0] != UINT_MAX && _fixed[k] != 1) {
      for (i = 0; i < movables[k].size(); ++i) {
        const HGFNode& node = hgwd.getNodeByIdx(movables[k][i]);
        for (itHGFEdgeLocal e = node.edgesBegin(); e != node.edgesEnd(); e++) {
          const HGFEdge& edge = (**e);
          unsigned netIdx = edge.getIndex();
          if (!_bitBoards.seenNets.isBitSet(netIdx)) {
            _bitBoards.seenNets.setBit(netIdx);
            unsigned haveFixed = false;
            BBox box;
            unsigned nodeOffset = 0;
            for (itHGFNodeLocal otherNode = edge.nodesBegin();
                 otherNode != edge.nodesEnd(); otherNode++, nodeOffset++) {
              unsigned nodeIndex = (*otherNode)->getIndex();
              if (_bitBoards.mobileCells.isBitSet(nodeIndex)) continue;

              if (placed[nodeIndex]) {
                //			     Point
                //pinOffset=hgwd.nodesOnEdgePinOffset
                //			       (nodeOffset,netIdx,pl.getOrient(nodeIndex));
                //			     box+=(pl[nodeIndex]+pinOffset);
                BBox pinOffset = hgwd.nodesOnEdgePinOffsetBBox(
                    nodeOffset, netIdx, pl.getOrient(nodeIndex));
                box += (pl[nodeIndex] + pinOffset.getBottomLeft());
                box += (pl[nodeIndex] + pinOffset.getTopRight());
              } else {
                box += pl[nodeIndex];
              }
              haveFixed = true;
            }
            // remove inessential nets
            if (!haveFixed || (box.xMin > rowsXMin || box.xMax < rowsXMax ||
                               box.yMin > rowsYMin || box.yMax < rowsYMax)) {
              _netTerminalBBoxes.push_back(box);
              _bitBoards.essentialNets.setBit(netIdx);
            }
          }
        }
      }
    }
  }
  const auto& essNetIds = _bitBoards.essentialNets.getIndicesOfSetBits();
  unsigned numNets = essNetIds.size();

  _netlist.setup(numCells, numNets);

  std::unordered_map<unsigned, unsigned, std::hash<unsigned>, std::equal_to<unsigned> >
      origToNewIdx;
  for (k = 0; k != numCells; k++)
    if (movables[k][0] != UINT_MAX)
      for (i = 0; i < movables[k].size(); ++i) origToNewIdx[movables[k][i]] = k;

  // setup the pin offsets for the movable cells
  for (k = 0; k != numNets; k++) {
    unsigned netIndex = essNetIds[k];
    const HGFEdge& edge = hgwd.getEdgeByIdx(netIndex);
    unsigned nodeOffset = 0;
    for (itHGFNodeLocal node = edge.nodesBegin(); node != edge.nodesEnd();
         node++, nodeOffset++) {
      unsigned nodeIndex = (*node)->getIndex();
      if (!_bitBoards.mobileCells.isBitSet(nodeIndex)) continue;

      // cout<<_bitBoards.clusteredCells.isBitSet(nodeIndex)<<endl;
      BBox pinOffset;
      if (_bitBoards.clusteredCells.isBitSet(nodeIndex) == 0) {
        //	      pinOffset += hgwd.nodesOnEdgePinOffset(nodeOffset,
        //	                                      netIndex,pl.getOrient(nodeIndex));
        BBox pinBBox = hgwd.nodesOnEdgePinOffsetBBox(nodeOffset, netIndex,
                                                     pl.getOrient(nodeIndex));
        pinOffset += pinBBox.getBottomLeft();
        pinOffset += pinBBox.getTopRight();
      } else  // clustered cell
      {
        bool found = 0;
        double leftEdge = 0;
        unsigned i = 0;
        for (i = 0; i < movables.size(); ++i) {
          found = 0;
          if (movables[i].size() == 1) continue;
          for (unsigned j = 0; j < movables[i].size(); ++j) {
            if (nodeIndex == movables[i][j]) {
              found = 1;
              break;
            }
          }
          if (found) break;
        }
        abkfatal(found, "clustered cell not found in list of cells\n");

        found = 0;
        leftEdge = 0;
        for (unsigned j = 0; j < movables[i].size(); ++j) {
          const HGFNode& otherNode = hgwd.getNodeByIdx(movables[i][j]);
          unsigned otherEdgeOffset = 0;
          for (itHGFEdgeLocal e = otherNode.edgesBegin();
               e != otherNode.edgesEnd(); e++, otherEdgeOffset++) {
            const HGFEdge& otherEdge = (**e);
            unsigned otherNetIdx = otherEdge.getIndex();
            if (otherNetIdx == netIndex) {
              //			  Point temp
              //=hgwd.edgesOnNodePinOffset(otherEdgeOffset, 			                netIndex,
              //pl.getOrient(movables[i][j]));
              //                          temp.x += leftEdge;
              //                          pinOffset += temp;
              BBox pinBBox = hgwd.edgesOnNodePinOffsetBBox(
                  otherEdgeOffset, netIndex, pl.getOrient(movables[i][j]));

              Point temp = pinBBox.getBottomLeft();
              temp.x += leftEdge;
              pinOffset += temp;

              temp = pinBBox.getTopRight();
              temp.x += leftEdge;
              pinOffset += temp;

              found = 1;
            }
          }
          leftEdge += hgwd.getNodeWidth(movables[i][j]);
        }
        abkfatal(found, "essential net not found in list of clustered cells\n")
      }
      _netlist.setPinOffset(origToNewIdx[nodeIndex], k, pinOffset);
    }
  }

  _bitBoards.seenNets.clear();
  _bitBoards.essentialNets.clear();
  _bitBoards.mobileCells.clear();
  _bitBoards.clusteredCells.clear();
}
