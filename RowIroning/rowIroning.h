/**************************************************************************
***    
*** Copyright (c) 2026 IonQ, Inc.
*** Copyright (c) 1995-2000 Regents of the University of California,
***               Andrew E. Caldwell, Andrew B. Kahng and Igor L. Markov
*** Copyright (c) 2000-2012 Regents of the University of Michigan,
***               Saurabh N. Adya, Jarrod A. Roy, David A. Papa and
***               Igor L. Markov
***
***  Contact author(s): abk@cs.ucsd.edu, imarkov@umich.edu
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


#ifndef  _ROW_IRONING_H_
#define  _ROW_IRONING_H_

#include "ABKCommon/abkcommon.h"
#include "HGraphWDims/hgWDims.h"
#include "RBPlace/RBPlacement.h"
#include "SmallPlacement/smallPlaceProb.h"
#include "Geoms/ptsetpt.h"
#include "riParams.h"

#include <iostream>

std::ostream& operator<<(std::ostream&, const RowIroningParameters&);

void ironRows(RBPlacement& rbplace, SmallPlacementBitBoardContainer &spBitBoards, MaxMem *maxMem);
void ironRows(RBPlacement& rbplace, const RowIroningParameters& params, SmallPlacementBitBoardContainer &spBitBoards, MaxMem *maxMem);

//added by DAP to allow parameter based calling of the LR and RL
//versions of this function
void ironRows2D(RBPlacement& rbplace, const RowIroningParameters& params, SmallPlacementBitBoardContainer &spBitBoards, MaxMem *maxMem,
				bool verticalDir, bool twoDim, bool cluster, bool leftToRight);

//added by sadya to form 2D instances for Rowironing
//allow left-right and right-left traversal
void ironRows2DRL(RBPlacement& rbplace, const RowIroningParameters& params, SmallPlacementBitBoardContainer &spBitBoards, MaxMem *maxMem,
				  bool verticalDir, bool twoDim, bool cluster);
void ironRows2DLR(RBPlacement& rbplace, const RowIroningParameters& params, SmallPlacementBitBoardContainer &spBitBoards, MaxMem *maxMem,
				  bool verticalDir, bool twoDim, bool cluster);

bool divideWhitespace(const RBPSubRow &subrow, unsigned windowSize,
                      bool LR, SmallPlacementRow &smallRow, double currWSWidth, Point currWSLoc,
                      unsigned &movablesCounter, std::vector< std::vector<unsigned> > &movables,
                      std::vector<double> &whiteSpaceWidths, std::vector<Point> &whiteSpaceLocs,
                      unsigned wsPieceCap);

#ifdef USEFLUTE
void placeWithFlute(const RBPlacement &rbplace, 
                    const uofm::vector< uofm::vector<unsigned> > &movables,
                    const SmallPlacementRow &row,
                    const uofm::vector<double> &whiteSpaceWidths,
                    BitBoard &nets, Placement &soln, Placement &currSoln,
                    uofm::vector< uofm::vector<PtSetPoint> > &ptSets);
void placeWithFluteDP(const RBPlacement &rbplace,
                      const uofm::vector< uofm::vector<unsigned> > &movables,
                      const SmallPlacementRow &row,
                      const uofm::vector<double> &whiteSpaceWidths,
                      BitBoard &nets, Placement &soln, Placement &tempSoln,
                      uofm::vector< uofm::vector<PtSetPoint> > &ptSets);

void addToPointSet(uofm::vector<PtSetPoint> &ptSet, const Point &pt);
void removeFromPointSet(uofm::vector<PtSetPoint> &ptSet, const Point &pt);
#endif

#endif
