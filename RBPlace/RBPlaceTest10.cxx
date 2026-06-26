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

// Placement-difference regression driver.
// It reads a saved placement, compares it to the current RBPlacement, and
// writes displacement vectors for cells that moved more than the threshold.

#include "RBPlacement.h"

using std::cout;
using std::endl;
using std::ifstream;
using std::ofstream;
using std::string;

int main(int argc, const char** argv) {
  StringParam auxFile("f", argc, argv);
  StringParam oldPlace("oldPlace", argc, argv);

  abkfatal(auxFile.found(), "-f auxFileName is required");
  abkfatal(oldPlace.found(), "-oldPlace placementName is required");

  RBPlacement::Parameters rbParams(argc, argv);
  RBPlacement rbplace(auxFile, rbParams);

  PlacementWOrient oldPl(rbplace.getHGraph().getNumNodes(),
                         Point(DBL_MAX, DBL_MAX));
  ifstream plFile(oldPlace);

  abkfatal2(plFile, " Could not open ", oldPlace);

  char nodeName[256];
  Point pt;
  char oriStr[10];
  int lineno = 0;

  plFile >> needcaseword("UCLA") >> needcaseword("pl") >> needword("1.0") >>
      skiptoeol;
  plFile >> eathash(lineno);

  while (!plFile.eof()) {
    plFile >> eathash(lineno) >> nodeName >> pt.x >> pt.y >> eatblank;
    unsigned nodeId = rbplace.getHGraph().getNodeByName(nodeName).getIndex();

    oldPl[nodeId] = pt;
    if (plFile.peek() == ':' || plFile.peek() == '/') {
      plFile.get();
      eatblank(plFile);
      if (plFile.peek() == '\n' || plFile.peek() == '\r') {
      } else {
        plFile >> oriStr;
        oldPl.setOrient(nodeId, Orientation(oriStr));
      }
      eatblank(plFile);
      if (plFile.peek() == ':' || plFile.peek() == '/') {
        plFile.get();
        eatblank(plFile);
        if (plFile.peek() == '\n' || plFile.peek() == '\r') {
        } else {
          string fixedStr;
          plFile >> fixedStr;
        }
      }
    }
    plFile >> skiptoeol >> eathash(lineno);
  }

  BBox core = rbplace.getBBox(true);
  double cutOffDist = 0.015 * (core.getWidth() + core.getHeight());
  unsigned displacedCells = 0;

  ofstream outFile("diffs.dat");

  double totalDist = 0.;
  double maxDist = -DBL_MAX;

  for (unsigned i = rbplace.getHGraph().getNumTerminals();
       i < rbplace.getHGraph().getNumNodes(); ++i) {
    Point newPlacement = rbplace[i];
    Point oldPlacement = oldPl[i];

    double dist = fabs(newPlacement.x - oldPlacement.x) +
                  fabs(newPlacement.y - oldPlacement.y);

    maxDist = std::max(maxDist, dist);
    totalDist += dist;

    if (greaterThanDouble(dist, cutOffDist)) {
      ++displacedCells;
      outFile << newPlacement << endl << oldPlacement << endl << endl;
    }
  }
  cout << "number of displacements plotted is " << displacedCells << endl;

  cout << "maxdist was " << maxDist << " which is "
       << 100. * maxDist / (core.getWidth() + core.getHeight())
       << "% of core HP" << endl;

  double avgDist =
      totalDist / static_cast<double>(rbplace.getHGraph().getNumNodes() -
                                      rbplace.getHGraph().getNumTerminals());

  cout << "average dist was " << avgDist << " which is "
       << 100. * avgDist / (core.getWidth() + core.getHeight())
       << "% of core HP" << endl;

  outFile.close();

  return 0;
}
