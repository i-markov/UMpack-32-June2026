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

// Placement width-adjustment driver.
// It rescales movable cell widths from the row-based benchmark and writes the
// updated design out for comparison.

#include "RBPlacement.h"

int main(int argc, const char **argv) {
  StringParam auxFile("f", argc, argv);
  StringParam save("saveAsNodes", argc, argv);
  DoubleParam setDims("setDims", argc, argv);
  StringParam newName("newName", argc, argv);

  abkfatal(auxFile.found(), "-f auxFileName is required");
  abkfatal(newName.found(), "-newName newFileName is required");

  RBPlacement::Parameters rbParams(argc, argv);
  RBPlacement rbplace(auxFile, rbParams);

  RandomRawUnsigned rng;

  double spacing = rbplace.getCoreRow(0).getSpacing();
  unsigned minWidth = UINT_MAX, maxWidth = 0;
  double rowHeight = rbplace.getCoreRow(0).getHeight();

  DoubleParam pctage("pctage", argc, argv);
  double cutoffPct = 0.;
  if (pctage.found()) {
    cutoffPct = 1. - pctage;
    cutoffPct = std::max(cutoffPct, 0.);
    cutoffPct = std::min(cutoffPct, 1.);
  }

  unsigned changeCutoff = static_cast<unsigned>(cutoffPct * (double)UINT_MAX);
  unsigned posNegCutoff = static_cast<unsigned>(0.5 * (double)UINT_MAX);

  unsigned total_increase = 0;
  unsigned total_decrease = 0;

  double areaBefore = rbplace.getMovableArea();

  HGraphWDimensions &hg = const_cast<HGraphWDimensions &>(rbplace.getHGraph());

  for (unsigned i = hg.getNumTerminals(); i < hg.getNumNodes(); ++i) {
    if (hg.getNodeHeight(i) <= rowHeight) {
      unsigned width = static_cast<unsigned>(hg.getNodeWidth(i) / spacing);
      minWidth = std::min(minWidth, width);
      maxWidth = std::max(maxWidth, width);
    }
  }

  for (unsigned i = hg.getNumTerminals(); i < hg.getNumNodes(); ++i) {
    if (hg.getNodeHeight(i) > rowHeight) {
      continue;
    }

    unsigned width = static_cast<unsigned>(hg.getNodeWidth(i) / spacing);
    unsigned range = std::min(width - minWidth, maxWidth - width);

    if (range > 0 && rng >= changeCutoff) {
      unsigned pos_neg = rng;
      unsigned magnitude = 1 + (rng % range);
      if (pos_neg > posNegCutoff) {
        width += magnitude;
        total_increase += magnitude;
      } else {
        width -= magnitude;
        total_decrease += magnitude;
      }
      hg.setNodeWidth(spacing * (double)width, i);
    }
  }

  //    cout << " number of cells " << hg.getNumNodes() - hg.getNumTerminals()
  //    << endl; cout << " area increase " << total_increase << endl; cout << "
  //    area decrease " << total_decrease << endl;

  rbplace.saveAsNodesNetsWts(newName);
  double areaAfter = rbplace.getMovableArea();
  double overlap = rbplace.calcOverlapGeneric(All);

  cout << endl;
  cout << "Area ratio is " << areaAfter / areaBefore << endl;
  cout << "Overlap% is " << 100. * overlap / areaAfter << endl;
  cout << endl;

  pair<double, double> HPWL = rbplace.evalXYHPWL(false);
  cout << " Center-to-center HalfPerim WL: " << HPWL.first + HPWL.second << " ("
       << HPWL.first << " x , " << HPWL.second << " y )" << endl;
  ;
  HPWL = rbplace.evalXYHPWL(true);
  cout << " Pin-to-pin HalfPerim WL      : " << HPWL.first + HPWL.second << " ("
       << HPWL.first << " x , " << HPWL.second << " y )" << endl;
  ;

  return 0;
}
