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

// created by Igor Markov on 06/02/98
#ifdef _MSC_VER
#pragma warning(disable : 4786)
#endif

#include "strayNodes.h"

#include "Partitioning/partitioning.h"
#include "partEvalHelpers.h"

using std::endl;
using std::ostream;
using std::setw;

StrayNodes::StrayNodes(const PartitioningProblem& problem,
                       const Partitioning& part, unsigned terminalsCountAs)
    : TalliesWCosts(problem, part, terminalsCountAs),
      _netDegree(problem.getHGraph().getNumEdges()),
      _maxPartSize(problem.getHGraph().getNumEdges()),
      _maxPart(problem.getHGraph().getNumEdges()),
      _maxNetCost(UINT_MAX) {
  reinitializeProper();
}

StrayNodes::StrayNodes(const PartitioningProblem& problem,
                       unsigned terminalsCountAs)
    : TalliesWCosts(problem, terminalsCountAs),
      _netDegree(problem.getHGraph().getNumEdges()),
      _maxPartSize(problem.getHGraph().getNumEdges()),
      _maxPart(problem.getHGraph().getNumEdges()),
      _maxNetCost(UINT_MAX) {}

void StrayNodes::reinitializeProper() {
  unsigned maxNetDegree = PartEvalHelpers::updateStrayNodesStats(
      _tallies, _nParts, _netDegree, _maxPartSize, _maxPart);
  _maxNetCost = maxNetDegree / 2;  // integer division
  updateAllCosts();
}

ostream& StrayNodes::prettyPrint(ostream& os) const {
  const UDenseMatrix& mat = _tallies;
  os << " NetId  Cost  Tally (term == " << terminalsCountAs()
     << " modules)  // netDegree  maxPartSize  maxPart " << endl;
  for (unsigned i = 0; i != mat.getNumRows(); i++) {
    os << setw(6) << i << ": " << setw(4) << getNetCost(i) << " ";
    for (unsigned j = 0; j != mat.getNumCols(); j++) {
      os << " " << setw(4) << mat(i, j);
    }
    os << "     // " << setw(4) << _netDegree[i] << setw(4) << _maxPartSize[i]
       << setw(4) << _maxPart[i] << endl;
  }
  os << "   Total cost : " << getTotalCost() << endl << endl;

  return os;
}
