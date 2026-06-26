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

// Regression test driver for Capo.
// It exercises representative package behavior with sample inputs or
// command-line options and prints results for the regression harness.
// Individual files vary the data set, but follow the same compare-against-
// baseline pattern.

// DB-backed placement driver that also accepts bin and hierarchy files.
// This variant is useful for testing Capo with a DB design plus explicit
// -bin and -hcl inputs, then comparing the reported wirelength metrics.
// Use -savePl to write the final placement.
// CHANGES

#ifdef _MSC_VER
#pragma warning(disable : 4786)
#endif

#include <iostream>

#include "ABKCommon/abkassert.h"
#include "ABKCommon/abkmessagebuf.h"
#include "ABKCommon/abktimer.h"
#include "ABKCommon/abkversion.h"
#include "RBPlace/RBPlacement.h"
#include "capoPlacer.h"
// #include "greedyOrientOpt.h"
// #include "rowdp.h"
// #include "rownet.h"
// #include "netopt.h"

int main(int argc, const char *argv[]) {
  NoParams noParams(argc, argv);
  BoolParam helpRequest("help", argc, argv);
  BoolParam helpRequest1("h", argc, argv);
  StringParam savePl("savePl", argc, argv);

  cout << "Capo 11.0 (June 2026)" << endl;
  cout << "Compiled " << __DATE__ << " at " << __TIME__ << endl << endl;

  cout
      << "\n(c) Copyright, 1998-2000 Regents of the University of California\n";
  cout << "(c) Copyright, 2000-2007 Regents of the University of Michigan\n\n";
  cout << " Authors: Saurabh N. Adya, Andrew E. Caldwell, Andrew B. Kahng\n";
  cout << "          Igor L. Markov, David A. Papa and Jarrod A. Roy\n";
  cout << " Email:   imarkov@eecs.umich.edu\n";
  cout << " http://vlsicad.eecs.umich.edu/BK/PDtools/" << endl << endl;

  cout << getABKMessageBuf() << flush;

  PRINT_VERSION_INFO

  cout << "Permission is hereby granted, without written agreement and \n"
       << "without license or royalty fee, to use, copy, modify, and \n"
       << "distribute and sell this software and its documentation for \n"
       << "any purpose, provided that the above copyright notice, this \n"
       << "permission notice, and the following two paragraphs appear \n"
       << "in all copies of this software as well as in all copies of \n"
       << "supporting documentation.\n\n";
  cout << "THIS SOFTWARE AND SUPPORTING DOCUMENTATION ARE PROVIDED \"AS IS\".\n"
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

  // ------------- Process parameters before getting bogged down in computations

  if (noParams.found() || helpRequest.found() || helpRequest1.found()) {
    cout << " -h		help" << endl;
    cout << " -help	help" << endl;
    cout << " -f      <auxFilename> auxfile for design to place" << endl;
    cout << " -bin  <binPlFileName> " << endl;
    cout << " -hcl    <hclFileName> " << endl;
    cout << " -savePl <fileName>	save resulting placement" << endl;

    return 0;
  }

  StringParam auxFileName("f", argc, argv);
  StringParam binFileName("bin", argc, argv);
  StringParam hclFileName("hcl", argc, argv);
  abkfatal(auxFileName.found(), "Usage: prog -f filename.aux <more params>");
  abkfatal(binFileName.found(), "Usage: prog -bin binfilename ");
  abkfatal(hclFileName.found(), "Usage: prog -hcl hclfilename ");

  RBPlaceParams rbParam(0, 0);
  RBPlacement rbplace(auxFileName, rbParam);

  bool saveMem = false;
  CapoParameters capoParams(argc, argv, saveMem);

  Timer capoTimer;
  CapoPlacer capo(rbplace, hclFileName, binFileName, capoParams);
  capoTimer.stop();

  cout << " Capo Runtime: " << capoTimer.getUserTime() << endl;

  cout << " Origin-to-Origin HPWL: " << rbplace.evalHPWL(false) << endl;
  ;
  cout << " Pin-to-Pin       WPWL: " << rbplace.evalHPWL(true) << endl;
  ;
  cout << " Origin-to-Origin  WWL: " << rbplace.evalWeightedWL(false) << endl;
  ;
  cout << " Pin-to-Pin        WWL: " << rbplace.evalWeightedWL(true) << endl;
  ;

  if (savePl.found()) rbplace.savePlacement(savePl);

  return 0;
}
