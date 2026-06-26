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

// Placement recovery utility used with Capo-generated .pl files.
// It reads a RowBasedPlacement benchmark, removes macro overlaps, marks the
// macros fixed, and writes out-recovered.pl for downstream placement checks.
// Use -f <filename.aux> to point it at a Bookshelf/RowBasedPlacement design.
#include <fstream>
#include <iostream>

#include "ABKCommon/pathDelims.h"
#include "Capo/capoPlacer.h"
#include "RBPlace/RBPlacement.h"
using std::cerr;
using std::cout;
using std::endl;

int main(int argc, const char *argv[]) {
  StringParam auxFileName("f", argc, argv);
  abkfatal(auxFileName.found(), "Usage: prog -f filename.aux \n");

  ifstream ifs(auxFileName);
  abkfatal3(ifs, "Cannot open ", auxFileName, "\n ");
  unsigned linesRead = 0;
  while (!ifs.eof()) {
    if (linesRead >= 10) {
      cerr << "\n No recognized format spec in first 10 lines of file\n    "
           << auxFileName << endl;
      cerr << " Currently supported formats are RowBasedPlacement" << endl;
      abkfatal(0, " ");
    }
    char buf[1000];
    ifs >> buf;
    ifs >> skiptoeol;
    linesRead++;
    if (buf[0] == '#') continue;
    buf[999] = '\0';
    if (strcasecmp(buf, "RowBasedPlacement") && strcasecmp(buf, "RBPlace"))
      abkfatal(0, "Only RowBasedPlacement is recognized by this program");

    // case 2: RBPlace
    cout << "Got RowBasedPlacement on input ..." << endl;
    RBPlaceParams rbParam(argc, argv);
    cout << rbParam << endl;
    RBPlacement rbplace(auxFileName, rbParam);

    // remove overlaps
    cout << "Removing overlaps" << endl;
    bool snapMacros = true;
    rbplace.removeMacroOverlaps(snapMacros);

    // fix macros to where they are
    cout << "Fixing macros" << endl;
    rbplace.markMacrosAsFixed();

    // write placement file
    cout << "Saving out-recovered.pl" << endl;
    bool labelFixed = true;
    rbplace.savePlacement("out-recovered.pl", labelFixed);
    break;
  }

  if (ifs.eof()) {
    cerr << "\n No recognized format spec in first 10 lines of file\n    "
         << auxFileName << endl;
    cerr << " Currently supported formats are RowBasedPlacement" << endl;
    abkfatal(0, " ");
  }
}
