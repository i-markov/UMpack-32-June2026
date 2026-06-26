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

// Row-based placement serialization driver.
// It loads an RBPlacement from -f and can write the design back out in
// multiple node/net/LEF-DEF formats for regression comparison.

#include <iostream>

#include "ABKCommon/abkassert.h"
#include "ABKCommon/paramproc.h"
#include "RBPlacement.h"

using std::cout;
using std::endl;

int main(int argc, const char *argv[]) {
  RBPlacement::Parameters rbParams(argc, argv);
  StringParam saveAsNodes("saveAsNodes", argc, argv);
  StringParam saveAsNodesShredHW("saveAsNodesShredHW", argc, argv);
  StringParam savePlacementUnShred("savePlacementUnShred", argc, argv);
  StringParam saveLEFDEF("saveLEFDEF", argc, argv);
  NoParams noParams(argc, argv);

  if (noParams.found()) {
    cout << " -f <auxfilename> \n"
         << " -saveLEFDEF  <basefilename> \n"
         << " -saveAsNodes <basefilename> " << endl
         << " -saveAsNodesShredHW <basefilename> " << endl
         << " -savePlacementUnShred  <filename> \n";
    return 0;
  }

  cout << "RBPlacement Parameters" << endl;
  cout << rbParams << endl;

  StringParam auxFile("f", argc, argv);

  abkfatal(auxFile.found(), "must have -f <auxfilename>");

  RBPlacement rbplace(auxFile, rbParams);

  if (saveAsNodes.found()) rbplace.saveAsNodesNetsWts(saveAsNodes);
  if (saveAsNodesShredHW.found())
    rbplace.saveAsNodesNetsWtsShredHW(saveAsNodesShredHW);
  if (savePlacementUnShred.found())
    rbplace.savePlacementUnShredHW(savePlacementUnShred);
  if (saveLEFDEF.found()) rbplace.saveLEFDEF(saveLEFDEF);

  return 0;
}
