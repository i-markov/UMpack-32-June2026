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

//   created by Igor Markov, April 30, 1999

#ifndef __SMALLPL_PROB_H__
#define __SMALLPL_PROB_H__

#include "Ctainers/bitBoard.h"
#include "Placement/placeWOri.h"
#include "baseSmallPlPr.h"

class SmallPlacementBitBoardContainer {
 public:
  BitBoard seenNets;
  BitBoard essentialNets;
  BitBoard mobileCells;
  // added by sadya to represent clustered cells
  BitBoard clusteredCells;

  SmallPlacementBitBoardContainer(void)
      : seenNets(200),
        essentialNets(200),
        mobileCells(50),
        clusteredCells(50) {}
};

class SmallPlacementProblem : public BaseSmallPlacementProblem {
  SmallPlacementBitBoardContainer& _bitBoards;

 protected:
  SmallPlacementProblem(unsigned numCells,
                        SmallPlacementBitBoardContainer& bitBoards,
                        Verbosity verb = Verbosity())
      : BaseSmallPlacementProblem(numCells, verb), _bitBoards(bitBoards) {}

 public:
  SmallPlacementProblem(const char* auxFile,
                        SmallPlacementBitBoardContainer& bitBoards,
                        Verbosity verb = Verbosity())
      : BaseSmallPlacementProblem(auxFile, verb), _bitBoards(bitBoards) {}

  SmallPlacementProblem(const HGraphWDimensions& hgraph,
                        const PlacementWOrient& placement,
                        const bit_vector& placed,
                        const std::vector<unsigned>& movables,
                        const SmallPlacementRow&,
                        SmallPlacementBitBoardContainer& bitBoards,
                        Verbosity verb = Verbosity("1 1 1"));

  // this constructor is added by saurabh adya to take care of fake cells(WS)
  // in the smallplacement problem
  SmallPlacementProblem(const HGraphWDimensions& hgraph,
                        const PlacementWOrient& placement,
                        const bit_vector& placed, const bit_vector& fixed,
                        const std::vector<unsigned>& movables,
                        const std::vector<SmallPlacementRow>&,
                        const std::vector<double>& whiteSpaceWidths,
                        SmallPlacementBitBoardContainer& bitBoards,
                        Verbosity verb = Verbosity("1 1 1"));

  // this constructor is added by saurabh adya to allow clustering of cells
  // in the smallplacement problem
  SmallPlacementProblem(const HGraphWDimensions& hgraph,
                        const PlacementWOrient& placement,
                        const bit_vector& placed, const bit_vector& fixed,
                        const std::vector<std::vector<unsigned> >& movables,
                        const std::vector<SmallPlacementRow>&,
                        const std::vector<double>& whiteSpaceWidths,
                        SmallPlacementBitBoardContainer& bitBoards,
                        Verbosity verb = Verbosity("1 1 1"));

  // for each node in hgraph, placement gives you the node's location,
  // orientations gives you its orientation, and placed tells you
  // if the cell has actually been placed or not.  If it has not
  // the nodes placed location will be used as the location of all
  // pins on the node.

  double getSnappedCellWidth(const HGraphWDimensions& hgwd, unsigned) const;
};

#endif
