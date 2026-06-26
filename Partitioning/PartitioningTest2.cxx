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

// Partitioning consistency driver.
// It loads the problem and checks that each node degree matches the number of
// incident edges reported by the hypergraph representation.

// CHANGES

#include "HGraph/hgFixed.h"
#include "ABKCommon/paramproc.h"
#include "partitioning.h"

using std::cout;
using std::endl;

int main(int argc, const char *argv[]) {
  BoolParam helpRequest1("h", argc, argv);
  BoolParam helpRequest2("help", argc, argv);

  StringParam baseFileName("f", argc, argv);
  PartitioningProblem::Parameters params(argc, argv);
  if (helpRequest1.found() || helpRequest2.found()) return 0;

  abkfatal(baseFileName.found(), " -f baseFileName not found");
  cout << params << endl;
  PartitioningProblem problem(baseFileName, params);

  for (itHGFNodeGlobal n = problem.getHGraph().nodesBegin();
       n != problem.getHGraph().nodesEnd(); n++) {
    unsigned checker = 0;
    for (itHGFEdgeLocal e = (*n)->edgesBegin(); e != (*n)->edgesEnd(); e++)
      checker++;

    if (checker != (*n)->getDegree()) {
      cout << " node " << (*n)->getIndex() << " degree " << (*n)->getDegree();
      cout << " actual " << checker << endl;
    }
  }

  return 0;
}
