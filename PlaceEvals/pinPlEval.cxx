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

#ifdef _MSC_VER
#pragma warning(disable : 4786)
#endif

#include <algorithm>
#include <climits>
#include <iostream>
#include <map>

#include "pinPlEval.h"

#include "ABKCommon/abkassert.h"
#include "ABKCommon/abktypes.h"
#include "ABKCommon/abkrand.h"
#ifdef USEFLUTE
#include "Flute/flute.h"
#endif
#include "GeomTrees/FastSteiner/global.h"
#include "GeomTrees/FastSteiner/random.h"
#include "GeomTrees/FastSteiner/stnr1.h"

using std::cout;
using std::endl;
using std::map;
using std::max;
using std::ostream;
using std::pair;
using std::string;
using std::vector;

namespace pinPlEval {

int set_of_b(unsigned idx, vector<int> tree) {
  int i = -90;
  while ((i = tree[idx]) != -1) idx = tree[idx];
  return idx;
}

void join_b(unsigned idx1, unsigned idx2, vector<int> &tree, bool coin) {
  int set1 = set_of_b(idx1, tree), set2 = set_of_b(idx2, tree);
  abkassert(tree[set1] == -1,
            "Invalid tree look up, Internal RBPlacePlotter Error");
  abkassert(tree[set2] == -1,
            "Invalid tree look up, Internal RBPlacePlotter Error");
  if (set1 != set2) {
    if (coin)
      tree[set1] = set2;
    else
      tree[set2] = set1;
  }
}

// the algorithm to find connected components in a graph
// TODO CITE SEDGEWICK ALGORITHMS PART 1-4, Ch. 1

double evalHPWL(const PlacementWOrient &pl, const HGraphWDimensions &hg,
                bool usePinLocs, bool centerOfGravityFSFPEval,
                bool pinToPinFSFPEval) {
  if (centerOfGravityFSFPEval)
    return evalFSFPCenterOfGravityHPWL(pl, hg, usePinLocs);
  else if (pinToPinFSFPEval)
    return evalFSFPPinToPinHPWL(pl, hg, usePinLocs);
  else
    return evalTraditionalHPWL(pl, hg, usePinLocs);
}

double evalTraditionalHPWL(const PlacementWOrient &pl,
                           const HGraphWDimensions &hg, bool usePinLocs) {
  double totalWL = 0.0;
  bool useWts = false;
  for (unsigned e = 0; e < hg.getNumEdges(); e++)
    totalWL += oneNetHPWL(pl, hg, e, usePinLocs, useWts);

  return totalWL;
}

double evalTraditionalWeightedHPWL(const PlacementWOrient &pl,
                                   const HGraphWDimensions &hg,
                                   bool usePinLocs) {
  double totalWL = 0.0;
  bool useWts = true;
  for (unsigned e = 0; e < hg.getNumEdges(); e++)
    totalWL += oneNetHPWL(pl, hg, e, usePinLocs, useWts);

  return totalWL;
}

// Function added by DAP, to be used in conjunction with
// Freeshape floorplanning.  It will compute the HPWL of nets
// not added during shredding with locations set to the
// actual pin locations of the nodes they connect
// Turn on with parameter -pinToPinFSFP
double evalFSFPPinToPinHPWL(const PlacementWOrient &pl,
                            const HGraphWDimensions &hg, bool usePinLocs) {
  double totalWL = 0.0;
  unsigned index = UINT_MAX;
  bool useWts = false;

  for (itHGFEdgeLocal e = hg.edgesBegin(); e != hg.edgesEnd(); e++) {
    index = (*e)->getIndex();

    // ignore all fake nets
    string netName(hg.getNetNameByIndex(index));
    if (netName.find("temp") != string::npos) continue;

    totalWL += oneNetHPWL(pl, hg, index, usePinLocs, useWts);
  }

  return totalWL;
}

// Function added by DAP, to be used in conjunction with
// Freeshape floorplanning.  It will compute the HPWL of nets
// not added during shredding with locations set to the center
// of gravity of the locations of the shreds
// Turn on with parameter -centerOfGravityFSFP
double evalFSFPCenterOfGravityHPWL(const PlacementWOrient &pl,
                                   const HGraphWDimensions &hg,
                                   bool usePinLocs) {
  double totalWL = 0.0;
  BBox netBBox;
  unsigned index;
  RandomRawUnsigned r;

  // compute connected components
  vector<int> tree(hg.getNumNodes(), -1);
  for (itHGFEdgeGlobal it = hg.edgesBegin(); it != hg.edgesEnd(); ++it) {
    string netName(hg.getNetNameByIndex((*it)->getIndex()));
    if (netName.find("temp") == string::npos) continue;
    const HGFEdge &this_net = *(*it);
    itHGFNodeLocal firstIt = this_net.nodesBegin();
    const HGFNode &firstNode = *(*firstIt);
    unsigned firstIdx = firstNode.getIndex();

    for (itHGFNodeLocal nodeIt = (*it)->nodesBegin();
         nodeIt != (*it)->nodesEnd(); ++nodeIt) {
      unsigned v = (*nodeIt)->getIndex();
      if (firstIdx != v) join_b(firstIdx, v, tree, r & 1u);
      // spam_tree(tree);
    }
  }

  // for each node,
  //  add its location to a list corresponding to its connected component
  map<unsigned, vector<Point> > component_points;
  for (unsigned i = 0; i < hg.getNumNodes(); i++)
    component_points[set_of_b(i, tree)].push_back(pl[i]);

  // for each connected component
  // compute centers of gravity and store them in a table
  map<unsigned, Point> center_of_gravity;
  for (map<unsigned, vector<Point> >::iterator it = component_points.begin();
       it != component_points.end(); ++it) {
    double N = (it->second).size();
    unsigned component_num = it->first;
    center_of_gravity[component_num] = Point(0, 0);
    for (vector<Point>::iterator it2 = (it->second).begin();
         it2 != (it->second).end(); ++it2) {
      center_of_gravity[component_num] += Point(it2->x / N, it2->y / N);
    }
  }
  // for each real net
  // for each node on the net
  // index the table appropriately and add its location to the bbox
  for (itHGFEdgeLocal e = hg.edgesBegin(); e != hg.edgesEnd(); e++) {
    netBBox.clear();
    index = (*e)->getIndex();

    // ignore all fake nets
    string netName(hg.getNetNameByIndex((*e)->getIndex()));
    if (netName.find("temp") != string::npos) continue;

    for (itHGFNodeLocal n = (*e)->nodesBegin(); n != (*e)->nodesEnd(); n++)
      netBBox += center_of_gravity[set_of_b((*n)->getIndex(), tree)];

    totalWL += netBBox.getHalfPerim();
  }

  return totalWL;
}

// by sadya to report X/Y components of HPWL
std::pair<double, double> evalXYHPWL(const PlacementWOrient &pl,
                                     const HGraphWDimensions &hg,
                                     bool usePinLocs) {
  double XWL = 0;
  double YWL = 0;
  std::pair<double, double> netWL;
  for (unsigned e = 0; e < hg.getNumEdges(); e++) {
    netWL = oneNetXYHPWL(pl, hg, e, usePinLocs);
    XWL += netWL.first;
    YWL += netWL.second;
  }
  return std::make_pair(XWL, YWL);
}

#define NumNetCostFactors 51

static double NetCostFactorTable[NumNetCostFactors] = {
    0.0,      1.000000, 1.000000, 1.000000, 1.082797, 1.153598, 1.220592,
    1.282322, 1.338501, 1.399094, 1.449260, 1.497380, 1.545500, 1.593619,
    1.641739, 1.689859, 1.730376, 1.770893, 1.811410, 1.851927, 1.892444,
    1.928814, 1.965184, 2.001553, 2.037923, 2.074293, 2.106117, 2.137941,
    2.169766, 2.201590, 2.233414, 2.264622, 2.295830, 2.327038, 2.358246,
    2.389454, 2.418676, 2.447898, 2.477119, 2.506341, 2.535563, 2.560954,
    2.586346, 2.611737, 2.637129, 2.662520, 2.688684, 2.714848, 2.741013,
    2.767177, 2.793341};

double evalWeightedWL(const PlacementWOrient &pl, const HGraphWDimensions &hg,
                      bool usePinLocs) {
  double totalWL = 0;
  bool useWts = false;
  for (unsigned e = 0; e < hg.getNumEdges(); e++) {
    double netsWL = oneNetHPWL(pl, hg, e, usePinLocs, useWts);
    unsigned edgeDegree = hg.getEdgeByIdx(e).getDegree();

    if (edgeDegree < NumNetCostFactors)
      netsWL *= NetCostFactorTable[edgeDegree];
    else
      netsWL *= NetCostFactorTable[NumNetCostFactors - 1];
    totalWL += netsWL;
  }

  return totalWL;
}

double evalNetsHPWL(const PlacementWOrient &pl, const HGraphWDimensions &hg,
                    const std::vector<unsigned> &netsToEval, bool usePinLocs) {
  double totalWL = 0;
  bool useWts = false;
  for (unsigned e = 0; e < netsToEval.size(); e++) {
    totalWL += oneNetHPWL(pl, hg, netsToEval[e], usePinLocs, useWts);
  }

  return totalWL;
}

// added by sadya to calculate the HPWL of a group of cells
//  will be used in overlap removal and rowIroning
double calcInstHPWL(const PlacementWOrient &pl, const HGraphWDimensions &hg,
                    std::vector<unsigned> cellIds) {
  bit_vector seenNets(hg.getNumEdges(), 0);

  double totalHPWL = 0;
  bool useWts = false;

  for (unsigned i = 0; i < cellIds.size(); ++i) {
    const HGFNode &node = hg.getNodeByIdx(cellIds[i]);
    for (itHGFEdgeLocal e = node.edgesBegin(); e != node.edgesEnd(); ++e) {
      const HGFEdge &edge = (**e);
      unsigned netIdx = edge.getIndex();
      if (!seenNets[netIdx]) {
        seenNets[netIdx] = 1;
        totalHPWL += oneNetHPWL(pl, hg, netIdx, true, useWts);
      }
    }
  }
  return totalHPWL;
}

// speed, by DP
double oneNetHPWL(const PlacementWOrient &pl, const HGraphWDimensions &hg,
                  unsigned netId, bool usePinLocs, bool useWts = false) {
  // this function has a lot of code repetition
  // because things have been manually inlined for speed concerns
  double HPWL = 0.0;

  vector<Pointf> pinLocs;

  // get pin locations
  const HGFEdge &net = hg.getEdgeByIdx(netId);
  itHGFNodeLocal nodeIt;
  unsigned nodeOffset;
  for (nodeIt = net.nodesBegin(), nodeOffset = 0; nodeIt != net.nodesEnd();
       nodeIt++, nodeOffset++) {
    const HGFNode &node = (**nodeIt);

    if (usePinLocs) {
      unsigned cellId = node.getIndex();
      BBox pinBBox = hg.nodesOnEdgePinOffsetBBox(nodeOffset, netId,
                                                 pl.getOrient(cellId));
      if (equalDouble(pinBBox.getWidth(), 0.) &&
          equalDouble(pinBBox.getHeight(), 0.)) {
        // only one pin
        Point pinLoc = pl[cellId];
        pinLoc.x += pinBBox.xMin;
        pinLoc.y += pinBBox.yMin;
        pinLocs.push_back(pinLoc);
      } else {
        // repeated pins
        Point pinLoc = pl[cellId];
        pinLoc.x += pinBBox.xMin;
        pinLoc.y += pinBBox.yMin;
        pinLocs.push_back(pinLoc);

        pinLoc = pl[cellId];
        pinLoc.x += pinBBox.xMax;
        pinLoc.y += pinBBox.yMax;
        pinLocs.push_back(pinLoc);
      }
    } else {
      unsigned cellId = node.getIndex();
      const unsigned &angle = pl.getOrient(cellId).getAngle();
      bool notRotateCell = angle == 0 || angle == 180;

      double cellHeight =
          (notRotateCell) ? hg.getNodeHeight(cellId) : hg.getNodeWidth(cellId);
      double cellWidth =
          (notRotateCell) ? hg.getNodeWidth(cellId) : hg.getNodeHeight(cellId);

      Point centerOffset(0.5 * cellWidth, 0.5 * cellHeight);
      Point pinLoc = pl[cellId] + centerOffset;
      pinLocs.push_back(pinLoc);
    }
  }

  unsigned nDegree = pinLocs.size();
  if (nDegree <= 1) {
    // one pin has 0 WL - do nothing
  } else if (nDegree == 2) {  // fast special case for 2 pin nets
    double pin0Locx = pinLocs[0].x;
    double pin0Locy = pinLocs[0].y;

    double pin1Locx = pinLocs[1].x;
    double pin1Locy = pinLocs[1].y;

    if (useWts)
      HPWL += net.getWeight() *
              (fabs(pin1Locx - pin0Locx) + fabs(pin1Locy - pin0Locy));
    else
      HPWL += fabs(pin1Locx - pin0Locx) + fabs(pin1Locy - pin0Locy);
  } else if (nDegree == 3) {
    double pin0Locx = pinLocs[0].x;
    double pin0Locy = pinLocs[0].y;

    double pin1Locx = pinLocs[1].x;
    double pin1Locy = pinLocs[1].y;

    double pin2Locx = pinLocs[2].x;
    double pin2Locy = pinLocs[2].y;

    double maxy = 0.0;
    double miny = 0.0;
    double maxx = 0.0;
    double minx = 0.0;
    if (pin1Locx > pin0Locx) {
      if (pin2Locx > pin1Locx) {
        maxx = pin2Locx;
        minx = pin0Locx;
      } else {
        maxx = pin1Locx;
        if (pin0Locx < pin2Locx)
          minx = pin0Locx;
        else
          minx = pin2Locx;
      }
    } else {
      if (pin2Locx > pin0Locx) {
        maxx = pin2Locx;
        minx = pin1Locx;
      } else {
        maxx = pin0Locx;
        if (pin1Locx < pin2Locx)
          minx = pin1Locx;
        else
          minx = pin2Locx;
      }
    }

    if (pin1Locy > pin0Locy) {
      if (pin2Locy > pin1Locy) {
        maxy = pin2Locy;
        miny = pin0Locy;
      } else {
        maxy = pin1Locy;
        if (pin0Locy < pin2Locy)
          miny = pin0Locy;
        else
          miny = pin2Locy;
      }
    } else {
      if (pin2Locy > pin0Locy) {
        maxy = pin2Locy;
        miny = pin1Locy;
      } else {
        maxy = pin0Locy;
        if (pin1Locy < pin2Locy)
          miny = pin1Locy;
        else
          miny = pin2Locy;
      }
    }

    if (useWts)
      HPWL += net.getWeight() * (maxx - minx + maxy - miny);
    else
      HPWL += maxx - minx + maxy - miny;
  } else  // degree >= 4
  {
    double maxx = 0.0;
    double maxy = 0.0;
    double minx = 0.0;
    double miny = 0.0;
    double halfPerim = 0.0;

    // Begin by initializing the locations to the first four pins
    // faster because it saves some comparing
    double pin0Locx = pinLocs[0].x;
    double pin0Locy = pinLocs[0].y;

    double pin1Locx = pinLocs[1].x;
    double pin1Locy = pinLocs[1].y;

    double pin2Locx = pinLocs[2].x;
    double pin2Locy = pinLocs[2].y;

    double pin3Locx = pinLocs[3].x;
    double pin3Locy = pinLocs[3].y;

    if (pin1Locx > pin0Locx) {
      if (pin3Locx > pin2Locx) {
        if (pin3Locx > pin1Locx)
          maxx = pin3Locx;
        else
          maxx = pin1Locx;

        if (pin2Locx > pin0Locx)
          minx = pin0Locx;
        else
          minx = pin2Locx;
      } else {
        if (pin2Locx > pin1Locx)
          maxx = pin2Locx;
        else
          maxx = pin1Locx;

        if (pin3Locx > pin0Locx)
          minx = pin0Locx;
        else
          minx = pin3Locx;
      }
    } else {
      if (pin3Locx > pin2Locx) {
        if (pin3Locx > pin0Locx)
          maxx = pin3Locx;
        else
          maxx = pin0Locx;

        if (pin2Locx > pin1Locx)
          minx = pin1Locx;
        else
          minx = pin2Locx;
      } else {
        if (pin2Locx > pin0Locx)
          maxx = pin2Locx;
        else
          maxx = pin0Locx;

        if (pin3Locx > pin1Locx)
          minx = pin1Locx;
        else
          minx = pin3Locx;
      }
    }

    if (pin1Locy > pin0Locy) {
      if (pin3Locy > pin2Locy) {
        if (pin3Locy > pin1Locy)
          maxy = pin3Locy;
        else
          maxy = pin1Locy;

        if (pin2Locy > pin0Locy)
          miny = pin0Locy;
        else
          miny = pin2Locy;
      } else {
        if (pin2Locy > pin1Locy)
          maxy = pin2Locy;
        else
          maxy = pin1Locy;

        if (pin3Locy > pin0Locy)
          miny = pin0Locy;
        else
          miny = pin3Locy;
      }
    } else {
      if (pin3Locy > pin2Locy) {
        if (pin3Locy > pin0Locy)
          maxy = pin3Locy;
        else
          maxy = pin0Locy;

        if (pin2Locy > pin1Locy)
          miny = pin1Locy;
        else
          miny = pin2Locy;
      } else {
        if (pin2Locy > pin0Locy)
          maxy = pin2Locy;
        else
          maxy = pin0Locy;

        if (pin3Locy > pin1Locy)
          miny = pin1Locy;
        else
          miny = pin3Locy;
      }
    }
    // end initialize

    // Loop over the rest of the pins in pairs (A,B)  starting at the fifth and
    // sixth pin compare nodes A and B, compare the greater to the max, and the
    // lesser to the min The loop will execute 0 times if there are exactly 5
    // pins
    unsigned A = 0;
    unsigned B = 0;
    unsigned numLoopIterations =
        (nDegree - 4) / 2;  // integer division floor((nDegree-4)/2.0);

    for (unsigned loopCt = 0; loopCt < numLoopIterations; ++loopCt) {
      A = loopCt * 2 + 0 + 4;
      B = loopCt * 2 + 1 + 4;

      double pinALocx = pinLocs[A].x;
      double pinALocy = pinLocs[A].y;

      double pinBLocx = pinLocs[B].x;
      double pinBLocy = pinLocs[B].y;

      if (pinBLocx > pinALocx) {
        if (pinBLocx > maxx) maxx = pinBLocx;
        if (pinALocx < minx) minx = pinALocx;
      } else {
        if (pinALocx > maxx) maxx = pinALocx;
        if (pinBLocx < minx) minx = pinBLocx;
      }

      if (pinBLocy > pinALocy) {
        if (pinBLocy > maxy) maxy = pinBLocy;
        if (pinALocy < miny) miny = pinALocy;
      } else {
        if (pinALocy > maxy) maxy = pinALocy;
        if (pinBLocy < miny) miny = pinBLocy;
      }
    }

    // finally, if nDegree is odd, pick up the last compare
    if (nDegree % 2 == 1) {
      // Begin Compute location of pin N
      unsigned N = nDegree - 1;

      double pinNLocx = pinLocs[N].x;
      double pinNLocy = pinLocs[N].y;

      if (pinNLocx > maxx)
        maxx = pinNLocx;
      else if (pinNLocx < minx)
        minx = pinNLocx;
      if (pinNLocy > maxy)
        maxy = pinNLocy;
      else if (pinNLocy < miny)
        miny = pinNLocy;
    }

    halfPerim = maxx - minx + maxy - miny;
    if (useWts)
      HPWL += net.getWeight() * halfPerim;
    else
      HPWL += halfPerim;
  }

  return HPWL;
}

// by sadya to get X/Y components of HPWL
std::pair<double, double> oneNetXYHPWL(const PlacementWOrient &pl,
                                       const HGraphWDimensions &hg,
                                       unsigned netId, bool usePinLocs) {
  abkassert(netId < hg.getNumEdges(), "netID out of range");

  BBox netBox, pinOffsetBBox;

  const HGFEdge &net = hg.getEdgeByIdx(netId);

  itHGFNodeLocal n;
  unsigned nodeOffset = 0;
  for (n = net.nodesBegin(); n != net.nodesEnd(); n++, nodeOffset++) {
    unsigned cellId = (*n)->getIndex();
    if (!usePinLocs) {
      const unsigned &angle = pl.getOrient(cellId).getAngle();
      bool notRotateCell = angle == 0 || angle == 180;

      double cellHeight =
          (notRotateCell) ? hg.getNodeHeight(cellId) : hg.getNodeWidth(cellId);
      double cellWidth =
          (notRotateCell) ? hg.getNodeWidth(cellId) : hg.getNodeHeight(cellId);

      Point centerOffset(0.5 * cellWidth, 0.5 * cellHeight);
      netBox += (pl[cellId] + centerOffset);
    } else {
      pinOffsetBBox =
          hg.nodesOnEdgePinOffsetBBox(nodeOffset, netId, pl.getOrient(cellId));
      pinOffsetBBox.TranslateBy(pl[cellId]);
      netBox.expandToInclude(pinOffsetBBox);
    }
  }

  return std::make_pair(netBox.getWidth(), netBox.getHeight());
}

// added by MRG to find mean of pins attached to nets of this cell
Point calcMeanLoc(const PlacementWOrient &pl, const HGraphWDimensions &hg,
                  unsigned cellId, bool usePinLocs) {
  unsigned pinCount = 0;
  Point meanLoc(0, 0);
  bit_vector seenNets(hg.getNumEdges(), 0);

  // for each net attached to the cell
  const HGFNode &node = hg.getNodeByIdx(cellId);
  for (itHGFEdgeLocal e = node.edgesBegin(); e != node.edgesEnd(); e++) {
    const HGFEdge &edge = (**e);
    unsigned netIdx = edge.getIndex();
    if (!seenNets[netIdx]) {
      seenNets[netIdx] = 1;
      // for each pin on the net
      unsigned nodeOffset = 0;
      for (itHGFNodeLocal n = edge.nodesBegin(); n != edge.nodesEnd();
           n++, nodeOffset++) {
        unsigned id = (*n)->getIndex();
        // skip myself
        if (id == cellId) continue;
        if (!usePinLocs) {
          meanLoc += pl[id];
        } else {
          Point pinOffset =
              hg.nodesOnEdgePinOffset(nodeOffset, netIdx, pl.getOrient(id));
          meanLoc += pl[id] + pinOffset;
        }
        pinCount++;
      }
    }
  }

  if (pinCount == 0) {
    meanLoc = pl[cellId];
  } else {
    meanLoc.x /= pinCount;
    meanLoc.y /= pinCount;
  }
  return (meanLoc);
}

#ifdef USEFLUTE
double evalHPWL_Flute(const PlacementWOrient &pl, const HGraphWDimensions &hg,
                      bool usePinLocs) {
  readLUT();
  double totalWL = 0.0;
  for (unsigned e = 0; e < hg.getNumEdges(); ++e)
    totalWL += oneNetHPWL_Flute(pl, hg, e, usePinLocs);

  return totalWL;
}

double evalHPWL_Flautist(const PlacementWOrient &pl,
                         const HGraphWDimensions &hg, const vector<BBox> &obs,
                         bool usePinLocs) {
  readLUT();
  vector<unsigned> relevantobs(obs.size(), 0);
  for (unsigned i = 0; i < obs.size(); ++i) {
    relevantobs[i] = i;
  }
  double totalWL = 0.0;
  for (unsigned e = 0; e < hg.getNumEdges(); ++e)
    totalWL += oneNetHPWL_Flautist(pl, hg, e, usePinLocs, obs, relevantobs);

  return totalWL;
}

double oneNetHPWL_Flute(const PlacementWOrient &pl, const HGraphWDimensions &hg,
                        unsigned netId, bool usePinLocs,
                        vector<pair<Point, Point> > *edges) {
  const HGFEdge &edge = hg.getEdgeByIdx(netId);
  vector<Point> points;
  unsigned nodeOffset = 0;
  for (itHGFNodeLocal nodeIt = edge.nodesBegin(); nodeIt != edge.nodesEnd();
       ++nodeIt, ++nodeOffset) {
    const HGFNode &node = (**nodeIt);
    BBox pinBBox = hg.nodesOnEdgePinOffsetBBox(nodeOffset, netId,
                                               pl.getOrient(node.getIndex()));

    if (usePinLocs) {
      unsigned cellId = node.getIndex();
      if (equalDouble(pinBBox.xMin, pinBBox.xMax) &&
          equalDouble(pinBBox.yMin, pinBBox.yMax)) {
        // only one pin
        points.push_back(pl[cellId] + pinBBox.getBottomLeft());
      } else {
        // repeated pins
        points.push_back(pl[cellId] + pinBBox.getBottomLeft());
        points.push_back(pl[cellId] + pinBBox.getTopRight());
      }
    } else {
      unsigned cellId = node.getIndex();
      const unsigned &angle = pl.getOrient(cellId).getAngle();
      bool notRotateCell = angle == 0 || angle == 180;

      double cellHeight =
          (notRotateCell) ? hg.getNodeHeight(cellId) : hg.getNodeWidth(cellId);
      double cellWidth =
          (notRotateCell) ? hg.getNodeWidth(cellId) : hg.getNodeHeight(cellId);

      Point centerOffset(0.5 * cellWidth, 0.5 * cellHeight);
      points.push_back(pl[cellId] + centerOffset);
    }
  }

  sort(points.begin(), points.end());
  vector<Point>::iterator new_end = unique(points.begin(), points.end());
  points.erase(new_end, points.end());

  double flutewl;
  if (points.size() <= 1) {
    flutewl = 0.;

    if (edges != NULL) {
      edges->clear();
    }
  } else if (points.size() == 2) {
    flutewl = fabs(points[0].x - points[1].x) + fabs(points[0].y - points[1].y);

    if (edges != NULL) {
      edges->clear();
      edges->push_back(std::make_pair(points[0], points[1]));
    }
  } else if (points.size() == 3) {
    double minx, maxx, miny, maxy;

    if (points[0].x < points[1].x) {
      minx = points[0].x;
      maxx = points[1].x;
    } else {
      minx = points[1].x;
      maxx = points[0].x;
    }
    if (points[2].x < minx) {
      minx = points[2].x;
    } else if (points[2].x > maxx) {
      maxx = points[2].x;
    }

    if (points[0].y < points[1].y) {
      miny = points[0].y;
      maxy = points[1].y;
    } else {
      miny = points[1].y;
      maxy = points[0].y;
    }
    if (points[2].y < miny) {
      miny = points[2].y;
    } else if (points[2].y > maxy) {
      maxy = points[2].y;
    }

    flutewl = (maxx - minx) + (maxy - miny);

    if (edges != NULL) {
      double x[3];
      double y[3];
      for (unsigned i = 0; i < points.size(); ++i) {
        x[i] = points[i].x;
        y[i] = points[i].y;
      }
      Tree flutetree = flute(points.size(), x, y, ACCURACY);

      edges->clear();
      int branchnum = 2 * flutetree.deg - 2;
      for (int i = 0; i < branchnum; ++i) {
        int n = flutetree.branch[i].n;
        abkfatal(n < branchnum, "problem with Flute tree");
        if (i == n) continue;
        Point first(flutetree.branch[i].x, flutetree.branch[i].y);
        Point second(flutetree.branch[n].x, flutetree.branch[n].y);

        if (first != second) {
          edges->push_back(std::make_pair(first, second));
        }
      }

      free(flutetree.branch);
    }
  } else if (points.size() > MAXD) {
    BBox net;
    for (unsigned i = 0; i < points.size(); ++i) {
      net += points[i];
    }
    flutewl = net.getHalfPerim();

    if (edges != NULL) {
      edges->clear();
    }
  } else {
    double *x = new double[points.size()];
    double *y = new double[points.size()];
    for (unsigned i = 0; i < points.size(); ++i) {
      x[i] = points[i].x;
      y[i] = points[i].y;
    }
    Tree flutetree = flute(points.size(), x, y, ACCURACY);
    flutewl = flutetree.length;

    delete[] x;
    delete[] y;

    if (edges != NULL) {
      edges->clear();
      int branchnum = 2 * flutetree.deg - 2;
      for (int i = 0; i < branchnum; ++i) {
        int n = flutetree.branch[i].n;
        abkfatal(n < branchnum, "problem with Flute tree");
        if (i == n) continue;
        Point first(flutetree.branch[i].x, flutetree.branch[i].y);
        Point second(flutetree.branch[n].x, flutetree.branch[n].y);

        if (first != second) {
          edges->push_back(std::make_pair(first, second));
        }
      }
    }
    free(flutetree.branch);
  }

  return flutewl;
}

double oneNetHPWL_Flautist(const PlacementWOrient &pl,
                           const HGraphWDimensions &hg, unsigned netId,
                           bool usePinLocs, const vector<BBox> &obs,
                           const vector<unsigned> &relevantobs,
                           vector<pair<Point, Point> > *edges) {
  const HGFEdge &edge = hg.getEdgeByIdx(netId);
  vector<Point> points;
  unsigned nodeOffset = 0;
  for (itHGFNodeLocal nodeIt = edge.nodesBegin(); nodeIt != edge.nodesEnd();
       ++nodeIt, ++nodeOffset) {
    const HGFNode &node = (**nodeIt);
    BBox pinBBox = hg.nodesOnEdgePinOffsetBBox(nodeOffset, netId,
                                               pl.getOrient(node.getIndex()));

    if (usePinLocs) {
      unsigned cellId = node.getIndex();
      if (equalDouble(pinBBox.xMin, pinBBox.xMax) &&
          equalDouble(pinBBox.yMin, pinBBox.yMax)) {
        // only one pin
        points.push_back(pl[cellId] + pinBBox.getBottomLeft());
      } else {
        // repeated pins
        points.push_back(pl[cellId] + pinBBox.getBottomLeft());
        points.push_back(pl[cellId] + pinBBox.getTopRight());
      }
    } else {
      unsigned cellId = node.getIndex();
      const unsigned &angle = pl.getOrient(cellId).getAngle();
      bool notRotateCell = angle == 0 || angle == 180;

      double cellHeight =
          (notRotateCell) ? hg.getNodeHeight(cellId) : hg.getNodeWidth(cellId);
      double cellWidth =
          (notRotateCell) ? hg.getNodeWidth(cellId) : hg.getNodeHeight(cellId);

      Point centerOffset(0.5 * cellWidth, 0.5 * cellHeight);
      points.push_back(pl[cellId] + centerOffset);
    }
  }

  sort(points.begin(), points.end());
  vector<Point>::iterator new_end = unique(points.begin(), points.end());
  points.erase(new_end, points.end());

  double flutewl = DBL_MAX;

  if (points.size() > MAXD) {
    BBox net;
    for (unsigned i = 0; i < points.size(); ++i) {
      net += points[i];
    }
    flutewl = net.getHalfPerim();

    if (edges != NULL) {
      edges->clear();
    }
  } else {
    double *x = new double[points.size()];
    double *y = new double[points.size()];
    for (unsigned i = 0; i < points.size(); ++i) {
      x[i] = points[i].x;
      y[i] = points[i].y;
    }
    unsigned legal = 0;
    Tree flutetree =
        flautist(points.size(), x, y, ACCURACY, obs, relevantobs, legal);
    flutewl = flutetree.length;

    delete[] x;
    delete[] y;

    if (edges != NULL) {
      edges->clear();
      int branchnum = 2 * flutetree.deg - 2;
      for (int i = 0; i < branchnum; ++i) {
        int n = flutetree.branch[i].n;
        abkfatal(n < branchnum, "problem with Flute tree");
        if (i == n) continue;
        Point first(flutetree.branch[i].x, flutetree.branch[i].y);
        Point second(flutetree.branch[n].x, flutetree.branch[n].y);

        if (first != second) {
          edges->push_back(std::make_pair(first, second));
        }
      }
    }
    free(flutetree.branch);
  }

  return flutewl;
}
#endif

void getBestSteinerTree(const PlacementWOrient &pl, const HGraphWDimensions &hg,
                        unsigned edgeId, bool usePinLocs, bool cleanup,
                        vector<pair<Point, Point> > &bestTree) {
#ifdef USEFLUTE
  readLUT();
  vector<pair<Point, Point> > fluteEdges;
#endif
  vector<pair<Point, Point> > fastStEdges;

  unsigned edgeDegree = hg.getEdgeByIdx(edgeId).getDegree();

  const unsigned scratchSize = max(1U, 8 * edgeDegree);
  fastSteiner::stnr1_package_init(8 * edgeDegree, 0);

  vector<fastSteiner::Point> pt(scratchSize);
  vector<long> parent(scratchSize);

#ifdef USEFLUTE
  double fluteLen = oneNetHPWL_Flute(pl, hg, edgeId, usePinLocs, &fluteEdges);
  double fastStLen = oneNetHPWL_FastSteiner(pl, hg, edgeId, usePinLocs,
                                            &pt[0], &parent[0], &fastStEdges);
#else
  oneNetHPWL_FastSteiner(pl, hg, edgeId, usePinLocs, &pt[0], &parent[0],
                         &fastStEdges);
#endif

  bestTree.clear();
#ifdef USEFLUTE
  if (fluteEdges.size() > 0 && fluteLen <= fastStLen) {
    bestTree = fluteEdges;
  } else
#endif
      if (fastStEdges.size() > 0) {
    bestTree = fastStEdges;
  }

  if (cleanup) fastSteiner::stnr1_package_done();
}

void printSteinerTrees(ostream &out, const PlacementWOrient &pl,
                       const HGraphWDimensions &hg, bool usePinLocs,
                       bool cleanup) {
#ifdef USEFLUTE
  readLUT();
  vector<pair<Point, Point> > fluteEdges;
#endif
  vector<pair<Point, Point> > fastStEdges;

  unsigned largestDegree = 0;
  for (unsigned e = 0; e < hg.getNumEdges(); ++e) {
    largestDegree = max(largestDegree, hg.getEdgeByIdx(e).getDegree());
  }

  if (largestDegree == 0) return;

  fastSteiner::stnr1_package_init(8 * largestDegree, 0);

  vector<fastSteiner::Point> pt(8 * largestDegree);
  vector<long> parent(8 * largestDegree);

  for (unsigned e = 0; e < hg.getNumEdges(); ++e) {
#ifdef USEFLUTE
    double fluteLen = oneNetHPWL_Flute(pl, hg, e, usePinLocs, &fluteEdges);
    double fastStLen =
        oneNetHPWL_FastSteiner(pl, hg, e, usePinLocs, &pt[0], &parent[0],
                               &fastStEdges);
#else
    oneNetHPWL_FastSteiner(pl, hg, e, usePinLocs, &pt[0], &parent[0],
                           &fastStEdges);
#endif

#ifdef USEFLUTE
    if (fluteEdges.size() > 0 && fluteLen <= fastStLen) {
      out << "Edge id: " << e;
      out << " Edge degree: " << hg.getEdgeByIdx(e).getDegree();
      out << " Tree size: " << fluteEdges.size() << " :";
      for (unsigned i = 0; i < fluteEdges.size(); ++i) {
        out << " ( ( " << fluteEdges[i].first.x << " " << fluteEdges[i].first.y
            << " ) ( " << fluteEdges[i].second.x << " "
            << fluteEdges[i].second.y << " ) )";
      }
      out << endl;
    } else
#endif
        if (fastStEdges.size() > 0) {
      out << "Edge id: " << e;
      out << " Edge degree: " << hg.getEdgeByIdx(e).getDegree();
      out << " Tree size: " << fastStEdges.size() << " :";
      for (unsigned i = 0; i < fastStEdges.size(); ++i) {
        out << " ( ( " << fastStEdges[i].first.x << " "
            << fastStEdges[i].first.y << " ) ( " << fastStEdges[i].second.x
            << " " << fastStEdges[i].second.y << " ) )";
      }
      out << endl;
    }
  }

  if (cleanup) fastSteiner::stnr1_package_done();
}

class LessFSPoint {
 public:
  LessFSPoint() {}
  bool operator()(const fastSteiner::Point &a, const fastSteiner::Point &b) {
    return (a.x < b.x) || ((a.x == b.x) && (a.y < b.y));
  }
};

class EqualFSPoint {
 public:
  EqualFSPoint() {}
  bool operator()(const fastSteiner::Point &a, const fastSteiner::Point &b) {
    return (a.x == b.x) && (a.y == b.y);
  }
};

// requires that necessary arrays be allocated before calling, so
// that they only need to be allocated once
// also assumes that the setup functions for FastSteiner have already
// been called
double oneNetHPWL_FastSteiner(const PlacementWOrient &pl,
                              const HGraphWDimensions &hg, unsigned netId,
                              bool usePinLocs, fastSteiner::Point *pt,
                              long *parent,
                              vector<pair<Point, Point> > *edges) {
  const HGFEdge &edge = hg.getEdgeByIdx(netId);
  long n_terminals = edge.getDegree();

  unsigned i = 0;
  unsigned nodeOffset = 0;
  for (itHGFNodeLocal nodeIt = edge.nodesBegin(); nodeIt != edge.nodesEnd();
       ++nodeIt, ++nodeOffset, ++i) {
    const HGFNode &node = (**nodeIt);
    BBox pinBBox = hg.nodesOnEdgePinOffsetBBox(nodeOffset, netId,
                                               pl.getOrient(node.getIndex()));

    if (usePinLocs) {
      unsigned cellId = node.getIndex();
      if (equalDouble(pinBBox.xMin, pinBBox.xMax) &&
          equalDouble(pinBBox.yMin, pinBBox.yMax)) {
        // only one pin
        Point pinLoc = pl[cellId];
        pinLoc.x += pinBBox.xMin;
        pinLoc.y += pinBBox.yMin;
        pt[i].x = pinLoc.x;
        pt[i].y = pinLoc.y;
      } else {
        // repeated pins
        Point pinLoc = pl[cellId];
        pinLoc.x += pinBBox.xMin;
        pinLoc.y += pinBBox.yMin;
        pt[i].x = pinLoc.x;
        pt[i].y = pinLoc.y;
        ++n_terminals;
        ++i;
        pinLoc = pl[cellId];
        pinLoc.x += pinBBox.xMax;
        pinLoc.y += pinBBox.yMax;
        pt[i].x = pinLoc.x;
        pt[i].y = pinLoc.y;
      }
    } else {
      unsigned cellId = node.getIndex();
      const unsigned &angle = pl.getOrient(cellId).getAngle();
      bool notRotateCell = angle == 0 || angle == 180;

      double cellHeight =
          (notRotateCell) ? hg.getNodeHeight(cellId) : hg.getNodeWidth(cellId);
      double cellWidth =
          (notRotateCell) ? hg.getNodeWidth(cellId) : hg.getNodeHeight(cellId);

      Point centerOffset(0.5 * cellWidth, 0.5 * cellHeight);
      Point pinLoc = pl[cellId] + centerOffset;
      pt[i].x = pinLoc.x;
      pt[i].y = pinLoc.y;
    }
  }

  std::sort(pt, pt + n_terminals, LessFSPoint());
  fastSteiner::Point *new_end =
      std::unique(pt, pt + n_terminals, EqualFSPoint());
  n_terminals = new_end - pt;

  long max_points = 4 * n_terminals;
  long n_rounds, n_phases;
  long n_points = n_terminals;
  long flags = 0;
  long max_rounds = DEFAULT_MAX_ROUNDS;
  long max_phases_per_round = DEFAULT_PHASES_PER_ROUND;
  long max_stnr_per_round = DEFAULT_STNR_PER_ROUND;
  long max_stnr_per_phase = DEFAULT_STNR_PER_PHASE;
  double cut_off = DEFAULT_CUT_OFF;

  if (n_terminals <= 1) {
    return 0.;

    if (edges != NULL) {
      edges->clear();
    }
  }
  if (n_terminals == 2) {
    return fabs(pt[0].x - pt[1].x) + fabs(pt[0].y - pt[1].y);

    if (edges != NULL) {
      edges->clear();
      edges->push_back(
          std::make_pair(Point(pt[0].x, pt[0].y), Point(pt[1].x, pt[1].y)));
    }
  }
  if (n_terminals == 3) {
    double minx, maxx, miny, maxy;

    if (pt[0].x < pt[1].x) {
      minx = pt[0].x;
      maxx = pt[1].x;
    } else {
      minx = pt[1].x;
      maxx = pt[0].x;
    }
    if (pt[2].x < minx) {
      minx = pt[2].x;
    } else if (pt[2].x > maxx) {
      maxx = pt[2].x;
    }

    if (pt[0].y < pt[1].y) {
      miny = pt[0].y;
      maxy = pt[1].y;
    } else {
      miny = pt[1].y;
      maxy = pt[0].y;
    }
    if (pt[2].y < miny) {
      miny = pt[2].y;
    } else if (pt[2].y > maxy) {
      maxy = pt[2].y;
    }

    if (edges != NULL) {
      edges->clear();
      // only 3 points, so 1 call should be sufficient
      fastSteiner::stnr1(&n_points, max_points, pt, parent, &n_rounds,
                         max_rounds, &n_phases, max_phases_per_round,
                         max_stnr_per_round, max_stnr_per_phase, cut_off,
                         flags);

      for (long i = 0; i < n_points; ++i) {
        if (parent[i] == i) continue;

        Point first(pt[i].x, pt[i].y);
        Point second(pt[parent[i]].x, pt[parent[i]].y);

        if (first != second) {
          edges->push_back(std::make_pair(first, second));
        }
      }
    }

    return (maxx - minx) + (maxy - miny);
  }

  double stnr_len1 = fastSteiner::stnr1(
      &n_points, max_points, pt, parent, &n_rounds, max_rounds, &n_phases,
      max_phases_per_round, max_stnr_per_round, max_stnr_per_phase, cut_off,
      flags);

  if (edges != NULL) {
    edges->clear();

    for (long i = 0; i < n_points; ++i) {
      if (parent[i] == i) continue;

      Point first(pt[i].x, pt[i].y);
      Point second(pt[parent[i]].x, pt[parent[i]].y);

      if (first != second) {
        edges->push_back(std::make_pair(first, second));
      }
    }
  }

  n_points = n_terminals;

  std::reverse(pt, pt + n_points);

  double stnr_len2 = fastSteiner::stnr1(
      &n_points, max_points, pt, parent, &n_rounds, max_rounds, &n_phases,
      max_phases_per_round, max_stnr_per_round, max_stnr_per_phase, cut_off,
      flags);

  if (edges != NULL && stnr_len2 < stnr_len1) {
    edges->clear();

    for (long i = 0; i < n_points; ++i) {
      if (parent[i] == i) continue;

      Point first(pt[i].x, pt[i].y);
      Point second(pt[parent[i]].x, pt[parent[i]].y);

      if (first != second) {
        edges->push_back(std::make_pair(first, second));
      }
    }
  }

  return std::min(stnr_len1, stnr_len2);
}

// allocates and deallocates its own arrays and calls the necessary
// setup functions for FastSteiner
// does not call the cleanup functions of Fast Steiner,
// so it is safe to use inside other environments that have
// already been using FastSteiner
double oneNetHPWL_FastSteiner(const PlacementWOrient &pl,
                              const HGraphWDimensions &hg, unsigned netId,
                              bool usePinLocs) {
  const HGFEdge &edge = hg.getEdgeByIdx(netId);
  long n_terminals = edge.getDegree();
  vector<fastSteiner::Point> pt(8 * n_terminals);

  unsigned i = 0;
  unsigned nodeOffset = 0;
  for (itHGFNodeLocal nodeIt = edge.nodesBegin(); nodeIt != edge.nodesEnd();
       ++nodeIt, ++nodeOffset, ++i) {
    const HGFNode &node = (**nodeIt);
    BBox pinBBox = hg.nodesOnEdgePinOffsetBBox(nodeOffset, netId,
                                               pl.getOrient(node.getIndex()));

    if (usePinLocs) {
      unsigned cellId = node.getIndex();
      if (equalDouble(pinBBox.xMin, pinBBox.xMax) &&
          equalDouble(pinBBox.yMin, pinBBox.yMax)) {
        // only one pin
        Point pinLoc = pl[cellId];
        pinLoc.x += pinBBox.xMin;
        pinLoc.y += pinBBox.yMin;
        pt[i].x = pinLoc.x;
        pt[i].y = pinLoc.y;
      } else {
        // repeated pins
        Point pinLoc = pl[cellId];
        pinLoc.x += pinBBox.xMin;
        pinLoc.y += pinBBox.yMin;
        pt[i].x = pinLoc.x;
        pt[i].y = pinLoc.y;
        ++n_terminals;
        ++i;
        pinLoc = pl[cellId];
        pinLoc.x += pinBBox.xMax;
        pinLoc.y += pinBBox.yMax;
        pt[i].x = pinLoc.x;
        pt[i].y = pinLoc.y;
      }
    } else {
      unsigned cellId = node.getIndex();
      const unsigned &angle = pl.getOrient(cellId).getAngle();
      bool notRotateCell = angle == 0 || angle == 180;

      double cellHeight =
          (notRotateCell) ? hg.getNodeHeight(cellId) : hg.getNodeWidth(cellId);
      double cellWidth =
          (notRotateCell) ? hg.getNodeWidth(cellId) : hg.getNodeHeight(cellId);

      Point centerOffset(0.5 * cellWidth, 0.5 * cellHeight);
      Point pinLoc = pl[cellId] + centerOffset;
      pt[i].x = pinLoc.x;
      pt[i].y = pinLoc.y;
    }
  }

  if (n_terminals <= 1) {
    return 0.;
  }
  if (n_terminals == 2) {
    double wl = fabs(pt[0].x - pt[1].x) + fabs(pt[0].y - pt[1].y);
    return wl;
  }
  if (n_terminals == 3) {
    double minx, maxx, miny, maxy;

    if (pt[0].x < pt[1].x) {
      minx = pt[0].x;
      maxx = pt[1].x;
    } else {
      minx = pt[1].x;
      maxx = pt[0].x;
    }
    if (pt[2].x < minx) {
      minx = pt[2].x;
    } else if (pt[2].x > maxx) {
      maxx = pt[2].x;
    }

    if (pt[0].y < pt[1].y) {
      miny = pt[0].y;
      maxy = pt[1].y;
    } else {
      miny = pt[1].y;
      maxy = pt[0].y;
    }
    if (pt[2].y < miny) {
      miny = pt[2].y;
    } else if (pt[2].y > maxy) {
      maxy = pt[2].y;
    }

    return (maxx - minx) + (maxy - miny);
  }

  std::sort(&pt[0], &pt[0] + n_terminals, LessFSPoint());
  fastSteiner::Point *new_end =
      std::unique(&pt[0], &pt[0] + n_terminals, EqualFSPoint());
  n_terminals = new_end - &pt[0];

  long max_points = 4 * n_terminals;
  vector<long> parent(max_points);
  long n_rounds, n_phases;
  long flags = 0;
  long max_rounds = DEFAULT_MAX_ROUNDS;
  long max_phases_per_round = DEFAULT_PHASES_PER_ROUND;
  long max_stnr_per_round = DEFAULT_STNR_PER_ROUND;
  long max_stnr_per_phase = DEFAULT_STNR_PER_PHASE;
  double cut_off = DEFAULT_CUT_OFF;

  fastSteiner::stnr1_package_init(max_points, 0);

  long n_points = n_terminals;

  double stnr_len1 = fastSteiner::stnr1(
      &n_points, max_points, &pt[0], &parent[0], &n_rounds, max_rounds,
      &n_phases, max_phases_per_round, max_stnr_per_round, max_stnr_per_phase,
      cut_off, flags);

  n_points = n_terminals;

  std::reverse(&pt[0], &pt[0] + n_points);

  double stnr_len2 = fastSteiner::stnr1(
      &n_points, max_points, &pt[0], &parent[0], &n_rounds, max_rounds,
      &n_phases, max_phases_per_round, max_stnr_per_round, max_stnr_per_phase,
      cut_off, flags);

  return std::min(stnr_len1, stnr_len2);
}

// calls the fast steiner "done" function, so not safe for use within Capo
double evalHPWL_FastSteiner(const PlacementWOrient &pl,
                            const HGraphWDimensions &hg, bool usePinLocs) {
  unsigned largestDegree = 0;
  for (unsigned e = 0; e < hg.getNumEdges(); ++e) {
    largestDegree = max(largestDegree, hg.getEdgeByIdx(e).getDegree());
  }

  if (largestDegree == 0) return 0.;

  fastSteiner::stnr1_package_init(8 * largestDegree, 0);

  vector<fastSteiner::Point> pt(8 * largestDegree);
  vector<long> parent(8 * largestDegree);

  double totalWL = 0.0;
  for (unsigned e = 0; e < hg.getNumEdges(); ++e) {
    totalWL += oneNetHPWL_FastSteiner(pl, hg, e, usePinLocs, &pt[0],
                                      &parent[0]);
  }

  fastSteiner::stnr1_package_done();

  return totalWL;
}

// calls the fast steiner "done" function, so not safe for use within Capo
double evalHPWL_SteinerCombined(const PlacementWOrient &pl,
                                const HGraphWDimensions &hg, bool usePinLocs) {
#ifdef USEFLUTE
  readLUT();
#endif
  unsigned largestDegree = 0;
  for (unsigned e = 0; e < hg.getNumEdges(); ++e) {
    largestDegree = max(largestDegree, hg.getEdgeByIdx(e).getDegree());
  }

  if (largestDegree == 0) return 0.;

  fastSteiner::stnr1_package_init(8 * largestDegree, 0);

  vector<fastSteiner::Point> pt(8 * largestDegree);
  vector<long> parent(8 * largestDegree);

  double totalWL = 0.0;
  for (unsigned e = 0; e < hg.getNumEdges(); ++e) {
    double netHPWL = DBL_MAX;

#ifdef USEFLUTE
    netHPWL = std::min(netHPWL, oneNetHPWL_Flute(pl, hg, e, usePinLocs));
#endif

    netHPWL = std::min(netHPWL,
                       oneNetHPWL_FastSteiner(pl, hg, e, usePinLocs, &pt[0],
                                              &parent[0]));

    totalWL += netHPWL;
  }

  fastSteiner::stnr1_package_done();

  return totalWL;
}

}  // namespace pinPlEval
