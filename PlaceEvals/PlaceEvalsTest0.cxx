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

// Placement-evaluator regression driver.
// It loads the placement, layout boxes, and hypergraph from disk, then
// computes the selected placement cost so the evaluator registry can be
// checked against the baseline.

#include <fstream>
#include <iostream>

#include "ABKCommon/abktimer.h"
#include "ABKCommon/paramproc.h"
#include "HGraph/hgFixed.h"
#include "Placement/layoutBBoxes.h"
#include "Placement/placement.h"
#include "expectedBBox.h"
#include "placeEvals.h"

using std::cout;
using std::endl;
using std::flush;
using std::ifstream;

int main(int argc, const char *argv[]) {
  StringParam netDName("net", argc, argv);
  StringParam areName("are", argc, argv);
  StringParam blockName("blk", argc, argv);
  StringParam placeName("pl", argc, argv);
  BoolParam help("help", argc, argv);
  BoolParam shorthelp("h", argc, argv);

  if (help || shorthelp) {
    cout << "Usage: placeEvalsTest0 -net <netD file> -are <are file>\n"
         << "             -blk <block file> -pl <placement file>\n"
         << "             [-h | -help]\n";
    return 0;
  }

  abkfatal(netDName.found(), "-net <netD file> not found");
  abkfatal(areName.found(), "-are <are file> not found");
  abkfatal(blockName.found(), "-blk <block file> not found");
  abkfatal(placeName.found(), "-pl <placement file> not found");

  LayoutBBoxes regions(blockName);

  cout << "  generate placement...\n" << flush;
  ifstream plInput(placeName);
  abkfatal(plInput, "error in reading the placement\n");
  Placement placement(plInput);
  plInput.close();

  cout << "  generate HGraph...\n" << flush;
  HGraphFixed hgraph(netDName, areName);

  PlaceEvalType pleval(argc, argv);
  cout << "  evaluate placement with " << pleval << "...\n" << flush;
  Timer tm;
  cout << "Cost: " << pleval.computeCost(hgraph, regions, placement) << endl;
  tm.stop();
  cout << tm << endl;

  return 0;
}
