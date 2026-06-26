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

// Weighted-hypergraph statistics and export demo.
// It loads an HGraphWDimensions, optionally rescales node dimensions, and
// prints summary statistics before saving in alternative formats.
// Use it to exercise the command-line and export paths.

#include "ABKCommon/paramproc.h"
#include "hgWDims.h"

using std::cout;
using std::endl;

int main(int argc, const char** argv) {
  StringParam auxFile("f", argc, argv);
  StringParam save("saveAsNodes", argc, argv);
  DoubleParam setDims("setDims", argc, argv);

  abkfatal(auxFile.found(), "-f auxFileName is required");

  HGraphParameters hgParams(argc, argv);
  hgParams.makeAllSrcSnk = false;

  cout << "HGraph Parameters" << endl << hgParams << endl;

  HGraphWDimensions hg(auxFile, NULL, NULL, hgParams);

  if (setDims.found()) hg.setNodeDims(setDims);

  if (save.found()) hg.saveAsNodesNetsWts(save);

  return 0;
}
