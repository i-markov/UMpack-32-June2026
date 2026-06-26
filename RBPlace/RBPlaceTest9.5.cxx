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

// Driven-length based width-adjustment driver.
// It estimates net-driven cell lengths, adjusts widths around the median
// demand, and saves the resulting benchmark.

#include "PlaceEvals/pinPlEval.h"
#include "RBPlacement.h"

using std::cout;
using std::endl;
using std::max;
using std::min;
using std::pair;
using std::vector;

int main(int argc, const char **argv) {
  StringParam auxFile("f", argc, argv);
  StringParam newName("newName", argc, argv);

  abkfatal(auxFile.found(), "-f auxFileName is required");
  abkfatal(newName.found(), "-newName newFileName is required");

  RBPlacement::Parameters rbParams(argc, argv);
  RBPlacement rbplace(auxFile, rbParams);

  double spacing = rbplace.getCoreRow(0).getSpacing();
  double rowHeight = rbplace.getCoreRow(0).getHeight();

  cout << "site spacing is " << spacing << endl;
  cout << "row height is " << rowHeight << endl;

  vector<double> netLengths(rbplace.getHGraph().getNumEdges(), 0.);
  vector<double> drivenLengths(rbplace.getHGraph().getNumNodes(), 0.);

  HGraphWDimensions &hg = const_cast<HGraphWDimensions &>(rbplace.getHGraph());

  bool usePinLocs = true;
  for (unsigned i = 0; i < rbplace.getHGraph().getNumEdges(); ++i) {
    netLengths[i] = pinPlEval::oneNetHPWL_FastSteiner(rbplace.getPlacement(),
                                                      hg, i, usePinLocs);
  }

  double totalDrivenLength = 0.;
  unsigned totalDrivingCells = 0;

  for (unsigned i = hg.getNumTerminals(); i < hg.getNumNodes(); ++i) {
    if (hg.getNodeHeight(i) > rowHeight) {
      continue;
    }

    const HGFNode &node = hg.getNodeByIdx(i);
    for (itHGFEdgeLocal e = node.edgesBegin(); e != node.edgesEnd(); ++e) {
      if ((*e)->isNodeSrc(&node)) {
        drivenLengths[i] += netLengths[(*e)->getIndex()];
      }
    }

    totalDrivenLength += drivenLengths[i];
    if (drivenLengths[i] > 0.) ++totalDrivingCells;
  }

  if (totalDrivingCells == 0) {
    cout << "no driving cells, something is wrong" << endl;
    return 0;
  }

  double avgDrivenLength =
      totalDrivenLength / static_cast<double>(totalDrivingCells);

  cout << "avg length " << avgDrivenLength << endl;

  double medianDrivenLength;

  vector<double> temp = drivenLengths;

  if (totalDrivingCells % 2) {  // odd number of useful elements
    unsigned dist = temp.size() - 1 - totalDrivingCells / 2;
    nth_element(temp.begin(), temp.begin() + dist, temp.end());

    medianDrivenLength = temp[dist];
  } else {  // even number of useful elements
    unsigned dist = temp.size() - 1 - totalDrivingCells / 2;
    nth_element(temp.begin(), temp.begin() + dist, temp.end());
    nth_element(temp.begin(), temp.begin() + dist + 1, temp.end());

    medianDrivenLength = 0.5 * (temp[dist] + temp[dist + 1]);
  }

  cout << "median length " << medianDrivenLength << endl;
  cout << "total driving cells " << totalDrivingCells << endl;

  double areaBefore = rbplace.getMovableArea();

  unsigned increases = 0;
  unsigned decreases = 0;
  unsigned equals = 0;

  for (unsigned i = hg.getNumTerminals(); i < hg.getNumNodes(); ++i) {
    if (hg.getNodeHeight(i) > rowHeight) {
      continue;
    }

    if (drivenLengths[i] == 0.) continue;

    double ratioToMedian =
        1. + 1. * ((drivenLengths[i] / medianDrivenLength) - 1.);

    ratioToMedian = max(ratioToMedian, 0.5);
    ratioToMedian = min(ratioToMedian, 1.5);

    double newWidth;
    if (ratioToMedian > 1.)
      newWidth = ceil(hg.getNodeWidth(i) * ratioToMedian / spacing) * spacing;
    else
      newWidth = floor(hg.getNodeWidth(i) * ratioToMedian / spacing) * spacing;

    if (newWidth > hg.getNodeWidth(i)) ++increases;
    if (newWidth < hg.getNodeWidth(i)) ++decreases;
    if (newWidth == hg.getNodeWidth(i)) ++equals;

    hg.setNodeWidth(newWidth, i);
  }

  rbplace.saveAsNodesNetsWts(newName);
  double areaAfter = rbplace.getMovableArea();
  double overlap = rbplace.calcOverlapGeneric(All);

  cout << endl;
  cout << "cells increased " << increases << endl;
  cout << "cells still equal " << equals << endl;
  cout << "cells decreased " << decreases << endl;
  cout << "Area ratio is " << areaAfter / areaBefore << endl;
  cout << "Overlap% is " << 100. * overlap / areaAfter << endl;
  cout << endl;

  pair<double, double> HPWL = rbplace.evalXYHPWL(false);
  cout << " Center-to-center HalfPerim WL: " << HPWL.first + HPWL.second << " ("
       << HPWL.first << " x , " << HPWL.second << " y )" << endl;
  HPWL = rbplace.evalXYHPWL(true);
  cout << " Pin-to-pin HalfPerim WL      : " << HPWL.first + HPWL.second << " ("
       << HPWL.first << " x , " << HPWL.second << " y )" << endl;

  return 0;
}
