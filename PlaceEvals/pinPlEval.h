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

#ifndef PINPLEVAL_H
#define PINPLEVAL_H

#include <map>
#include <vector>

#include <string>
#include <vector>
#include "GeomTrees/FastSteiner/global.h"
#include "HGraphWDims/hgWDims.h"
#include "Placement/placeWOri.h"

namespace pinPlEval {
double evalHPWL(const PlacementWOrient &pl, const HGraphWDimensions &hg,
                bool usePinLocs, bool centerOfGravityFSFPEval,
                bool pinToPinFSFPEval);

double evalTraditionalHPWL(const PlacementWOrient &pl,
                           const HGraphWDimensions &hg, bool usePinLocs);

double evalTraditionalWeightedHPWL(const PlacementWOrient &pl,
                                   const HGraphWDimensions &hg,
                                   bool usePinLocs);

double evalFSFPPinToPinHPWL(const PlacementWOrient &pl,
                            const HGraphWDimensions &hg, bool usePinLocs);

double evalFSFPCenterOfGravityHPWL(const PlacementWOrient &pl,
                                   const HGraphWDimensions &hg,
                                   bool usePinLocs);

#ifdef USEFLUTE
double evalHPWL_Flute(const PlacementWOrient &pl, const HGraphWDimensions &hg,
                      bool usePinLocs);

double evalHPWL_Flautist(const PlacementWOrient &pl,
                         const HGraphWDimensions &hg,
                         const std::vector<BBox> &obs, bool usePinLocs);
#endif

double evalHPWL_FastSteiner(const PlacementWOrient &pl,
                            const HGraphWDimensions &hg, bool usePinLocs);

double evalHPWL_SteinerCombined(const PlacementWOrient &pl,
                                const HGraphWDimensions &hg, bool usePinLocs);

double evalWeightedWL(const PlacementWOrient &pl, const HGraphWDimensions &hg,
                      bool usePinLocs);

std::pair<double, double> evalXYHPWL(const PlacementWOrient &pl,
                                     const HGraphWDimensions &hg,
                                     bool usePinLocs);

double evalNetsHPWL(const PlacementWOrient &pl, const HGraphWDimensions &hg,
                    const std::vector<unsigned> &netsToEval, bool usePinLocs);

double oneNetHPWL(const PlacementWOrient &pl, const HGraphWDimensions &hg,
                  unsigned netId, bool usePinLocs, bool useWts);

#ifdef USEFLUTE
double oneNetHPWL_Flute(const PlacementWOrient &pl, const HGraphWDimensions &hg,
                        unsigned netId, bool usePinLocs,
                        std::vector<std::pair<Point, Point> > *edges = NULL);

double oneNetHPWL_Flautist(
    const PlacementWOrient &pl, const HGraphWDimensions &hg, unsigned netId,
    bool usePinLocs, const std::vector<BBox> &obs,
    const std::vector<unsigned> &relevantobs,
    std::vector<std::pair<Point, Point> > *edges = NULL);
#endif

double oneNetHPWL_FastSteiner(
    const PlacementWOrient &pl, const HGraphWDimensions &hg, unsigned netId,
    bool usePinLocs, fastSteiner::Point *pt, long *parent,
    std::vector<std::pair<Point, Point> > *edges = NULL);

double oneNetHPWL_FastSteiner(const PlacementWOrient &pl,
                              const HGraphWDimensions &hg, unsigned netId,
                              bool usePinLocs);

void getBestSteinerTree(const PlacementWOrient &pl, const HGraphWDimensions &hg,
                        unsigned edgeId, bool usePinLocs, bool cleanup,
                        std::vector<std::pair<Point, Point> > &bestTree);

void printSteinerTrees(std::ostream &out, const PlacementWOrient &pl,
                       const HGraphWDimensions &hg, bool usePinLocs,
                       bool cleanup);

std::pair<double, double> oneNetXYHPWL(const PlacementWOrient &pl,
                                       const HGraphWDimensions &hg,
                                       unsigned netId, bool usePinLocs);

double calcInstHPWL(const PlacementWOrient &pl, const HGraphWDimensions &hg,
                    std::vector<unsigned> cellIds);

Point calcMeanLoc(const PlacementWOrient &pl, const HGraphWDimensions &hg,
                  unsigned cellId, bool usePinLocs);
}  // namespace pinPlEval

#endif
