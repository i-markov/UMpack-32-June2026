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

#include "ABKCommon/abkassert.h"
#include "ABKCommon/abkversion.h"
#include "ABKCommon/infolines.h"
#include "ABKCommon/paramproc.h"
#include "ABKCommon/abkrand.h"
#include "PartEvals/partEvals.h"
#include "fmPart.h"

#include <iostream>

using std::cout;
using std::endl;

// Netlist export utility for weight-normalized inputs.
// Rewrites every hyperedge weight to a fixed value, reports the average
// weight used, and saves the modified design as a nodes/nets output.
int main(int argc, const char* argv[]) {
  BoolParam helpRequest1("h", argc, argv);
  BoolParam helpRequest2("help", argc, argv);

  StringParam baseFileName("f", argc, argv);
  UnsignedParam numStarts("num", argc, argv);

  BoolParam saveSoln("save", argc, argv);

  PRINT_VERSION_INFO

  PartitioningProblem::Parameters params(argc, argv);
  PartitionerParams pParams(argc, argv);
  if (helpRequest1.found() || helpRequest2.found()) return 0;
  abkfatal(baseFileName.found(), " -f baseFileName not found");
  cout << params << endl;
  PartitioningProblem problem(baseFileName, params);
  problem.discardHGraphNames();

  // change weights here

  RandomDouble r(80, 120);
  RandomUnsigned ii(0, problem.getHGraph().getNumEdges());

  double sum = 0;
  for (unsigned e = 0; e < problem.getHGraph().getNumEdges(); e++) {
    // double rr=r;
    double rr = 1000;
    sum += rint(rr);
    HGFEdge& ee = const_cast<HGFEdge&>(problem.getHGraph().getEdgeByIdx(e));
    ee.setWeight(rr);
  }
  cout << "The final average was" << (sum / problem.getHGraph().getNumEdges())
       << endl;

  // end weights change

  problem.saveAsNodesNets("newOUT");
}
