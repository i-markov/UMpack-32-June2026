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

// Layout-generation regression driver.
// It reads a benchmark and builds a new row-based layout using the requested
// aspect ratio and whitespace, then optionally saves the result.

#include <iostream>

#include "ABKCommon/paramproc.h"
#include "RBPlacement.h"

using std::cout;
using std::endl;

int main(int argc, const char *argv[]) {
  cout << "UCLA Layout Generator version 1.0  99/09/14" << endl;
  cout
      << "\n(c) copyright, 1998-1999 Regents of the University of California\n";
  cout << " Authors: Andrew E. Caldwell\n";
  cout << " Email:   caldwell@cs.ucla.edu" << endl << endl;
  cout << "Permission is hereby granted, without written agreement and \n"
       << "without license or royalty fee, to use, copy, modify, and \n"
       << "distribute and sell this software and its documentation for \n"
       << "any purpose, provided that the above copyright notice, this \n"
       << "permission notice, and the following two paragraphs appear \n"
       << "in all copies of this software as well as in all copies of \n"
       << "supporting documentation.\n\n";
  cout << "THIS SOFTWARE AND SUPPORTING DOCUMENTATION IS PROVIDED \"AS IS\".\n"
       << "The Microelectronics Advanced Research Corporation (MARCO), the\n"
       << "Gigascale Silicon Research Center (GSRC),  and \n"
       << "(\"PROVIDERS\") MAKE NO WARRANTIES, whether express \n"
       << "or implied, including warranties of merchantability or fitness\n"
       << "for a particular purpose or noninfringement, with respect to \n"
       << "this software and supporting documentation.\n"
       << "Providers have NO obligation to provide ANY support, assistance,\n"
       << "installation, training or other services, updates, enhancements\n"
       << "or modifications related to this software and supporting\n"
       << "documentation.\n\n";
  cout << "Providers shall NOT be liable for ANY costs of procurement of \n"
       << "substitutes, loss of profits, interruption of business, or any \n"
       << "other direct, indirect, special, consequential or incidental \n"
       << "damages arising from the use of this software and its \n"
       << "documentation, whether or not Providers have been advised of \n"
       << "the possibility of such damages." << endl
       << endl;

  BoolParam help1("h", argc, argv);
  BoolParam help2("help", argc, argv);
  StringParam auxFile("f", argc, argv);

  if (help1.found() || help2.found() || !auxFile.found()) {
    cout << "  Generates layouts for netlist in netD/areM or \n"
         << "  nodes/nets/wts formats." << endl;
    cout << "    -f <auxFileName>" << endl;
    cout << "    -AR <aspect ratio> (default is 1.0)" << endl;
    cout << "    -WS <whitespace %> (default is 10)" << endl;
    cout << "    -saveAs <newFileName> " << endl;

    return 0;
  }

  RBPlacement::Parameters rbParams(argc, argv);
  DoubleParam AR("AR", argc, argv);
  DoubleParam WS("WS", argc, argv);
  StringParam saveAs("saveAs", argc, argv);

  double ws;
  if (WS.found())
    ws = WS;
  else
    ws = 10;

  double ar;
  if (AR.found()) {
    ar = AR;
    cout << "Setting Aspect Ratio to " << ar << endl;
  } else {
    ar = 1;
    cout << "Using default Aspect Ratio of " << ar << endl;
  }

  RBPlacement rbplace(auxFile, ar, ws, rbParams);

  if (saveAs.found()) rbplace.saveAsNodesNetsWts(saveAs);

  return 0;
}
