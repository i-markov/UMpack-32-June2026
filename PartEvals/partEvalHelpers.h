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

#ifndef __PART_EVAL_HELPERS_H__
#define __PART_EVAL_HELPERS_H__

#include <iomanip>
#include <ostream>
#include <vector>

#include "Ctainers/umatrix.h"
#include "HGraph/hgFixed.h"

namespace PartEvalHelpers {

inline std::ostream& printTalliesWCosts2way(
    std::ostream& os, const HGraphFixed& hg, const unsigned short* tallies,
    unsigned terminalsCountAs, const std::vector<unsigned>& netCosts,
    const std::vector<unsigned>* netWeights, unsigned totalCost) {
  if (netWeights)
    os << "  NetId    Cost   Weight      Tally (terminals count as "
       << terminalsCountAs << " module(s)) " << std::endl;
  else
    os << "  NetId    Cost       Tally (terminals count as "
       << terminalsCountAs << " module(s)) " << std::endl;

  for (unsigned i = 0; i != hg.getNumEdges(); i++) {
    os << std::setw(6) << i << " : " << std::setw(6) << netCosts[i] << " :  ";
    if (netWeights) os << std::setw(3) << (*netWeights)[i] << " :  ";
    os << " " << std::setw(4) << static_cast<unsigned>(tallies[2 * i]) << " "
       << std::setw(4) << static_cast<unsigned>(tallies[2 * i + 1])
       << std::endl;
  }
  os << "   Total cost : " << totalCost << std::endl << std::endl;

  return os;
}

inline unsigned updateStrayNodesStats(const UDenseMatrix& tallies,
                                      unsigned nParts,
                                      std::vector<unsigned>& netDegree,
                                      std::vector<unsigned>& maxPartSize,
                                      std::vector<unsigned>& maxPart) {
  unsigned maxNetDegree = 0;
  for (unsigned netIndex = 0; netIndex != tallies.getNumRows(); netIndex++) {
    unsigned curMaxPartSize = tallies(netIndex, 0);
    unsigned curMaxPart = 0;
    unsigned curDegree = curMaxPartSize;

    for (unsigned part = 1; part != nParts; part++) {
      unsigned partSize = tallies(netIndex, part);
      if (partSize > curMaxPartSize) {
        curMaxPartSize = partSize;
        curMaxPart = part;
      }
      curDegree += partSize;
    }
    netDegree[netIndex] = curDegree;
    if (curDegree > maxNetDegree) maxNetDegree = curDegree;
    maxPartSize[netIndex] = curMaxPartSize;
    maxPart[netIndex] = curMaxPart;
  }
  return maxNetDegree;
}

inline unsigned maxNetDegree2way(const HGraphFixed& hg,
                                 const unsigned short* tallies) {
  unsigned maxNetDegree = 0;
  for (const unsigned short* pt = tallies; pt != tallies + 2 * hg.getNumEdges();) {
    unsigned deg = *pt++;
    deg += *pt++;
    if (deg > maxNetDegree) maxNetDegree = deg;
  }
  return maxNetDegree;
}

}  // namespace PartEvalHelpers

#endif
