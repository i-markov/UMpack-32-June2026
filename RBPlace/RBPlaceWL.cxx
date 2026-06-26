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

#ifdef _MSC_VER
#pragma warning(disable : 4786)
#endif

#include <iostream>
#include <map>

#include "PlaceEvals/pinPlEval.h"
#include "RBPlacement.h"

using std::cout;
using std::endl;
using std::map;
using std::string;

#ifdef _MSC_VER
#ifndef rint
#define rint(a) floor((a) + 0.5)
#endif
#endif

double RBPlacement::evalHPWL(const PlacementWOrient &pl,
                             const HGraphWDimensions &hg,
                             bool usePinLocs) const {
  return pinPlEval::evalHPWL(pl, hg, usePinLocs,
                             _params.centerOfGravityFSFPEval,
                             _params.pinToPinFSFPEval);
}

double RBPlacement::evalTraditionalHPWL(const PlacementWOrient &pl,
                                        const HGraphWDimensions &hg,
                                        bool usePinLocs) const {
  return pinPlEval::evalTraditionalHPWL(pl, hg, usePinLocs);
}

double RBPlacement::evalTraditionalWeightedHPWL(const PlacementWOrient &pl,
                                                const HGraphWDimensions &hg,
                                                bool usePinLocs) const {
  return pinPlEval::evalTraditionalWeightedHPWL(pl, hg, usePinLocs);
}

double RBPlacement::evalFSFPPinToPinHPWL(const PlacementWOrient &pl,
                                         const HGraphWDimensions &hg,
                                         bool usePinLocs) const {
  return pinPlEval::evalFSFPPinToPinHPWL(pl, hg, usePinLocs);
}

double RBPlacement::evalFSFPCenterOfGravityHPWL(const PlacementWOrient &pl,
                                                const HGraphWDimensions &hg,
                                                bool usePinLocs) const {
  return pinPlEval::evalFSFPCenterOfGravityHPWL(pl, hg, usePinLocs);
}

std::pair<double, double> RBPlacement::evalXYHPWL(const PlacementWOrient &pl,
                                                  const HGraphWDimensions &hg,
                                                  bool usePinLocs) const {
  return pinPlEval::evalXYHPWL(pl, hg, usePinLocs);
}

double RBPlacement::evalWeightedWL(const PlacementWOrient &pl,
                                   const HGraphWDimensions &hg,
                                   bool usePinLocs) const {
  return pinPlEval::evalWeightedWL(pl, hg, usePinLocs);
}

double RBPlacement::evalNetsHPWL(const PlacementWOrient &pl,
                                 const HGraphWDimensions &hg,
                                 const std::vector<unsigned> &netsToEval,
                                 bool usePinLocs) const {
  return pinPlEval::evalNetsHPWL(pl, hg, netsToEval, usePinLocs);
}

double RBPlacement::oneNetHPWL(const PlacementWOrient &pl,
                               const HGraphWDimensions &hg, unsigned netId,
                               bool usePinLocs) const {
  bool useWts = false;
  return pinPlEval::oneNetHPWL(pl, hg, netId, usePinLocs, useWts);
}

std::pair<double, double> RBPlacement::oneNetXYHPWL(const PlacementWOrient &pl,
                                                    const HGraphWDimensions &hg,
                                                    unsigned netId,
                                                    bool usePinLocs) {
  return pinPlEval::oneNetXYHPWL(pl, hg, netId, usePinLocs);
}

double RBPlacement::calcInstHPWL(const PlacementWOrient &pl,
                                 const HGraphWDimensions &hg,
                                 std::vector<unsigned> cellIds) {
  return pinPlEval::calcInstHPWL(pl, hg, cellIds);
}

Point RBPlacement::calcMeanLoc(const PlacementWOrient &pl,
                               const HGraphWDimensions &hg, unsigned cellId,
                               bool usePinLocs) {
  return pinPlEval::calcMeanLoc(pl, hg, cellId, usePinLocs);
}
