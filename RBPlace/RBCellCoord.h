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

#ifndef _RBCellCoord_H_
#define _RBCellCoord_H_

#include <climits>
#include <ostream>

// A coordinate 3-tuple identifies a cell by
//   2: RBPlacement coreRow index
//   3: RBPCoreRow subRow index
//   4: RBPSubRow cell index
// This allows finding a cell's location within a subRow in O(1) time

class RBCellCoord {
 public:
  unsigned coreRowId;   // RBPlacement coreRow index
  unsigned subRowId;    // RBPCoreRow subRow index
  unsigned cellOffset;  // Cell's position within the subRow

  RBCellCoord()
      : coreRowId(UINT_MAX), subRowId(UINT_MAX), cellOffset(UINT_MAX) {}

  RBCellCoord(unsigned crId, unsigned srId, unsigned cOffset)
      : coreRowId(crId), subRowId(srId), cellOffset(cOffset) {}

  RBCellCoord(const RBCellCoord& orig)
      : coreRowId(orig.coreRowId),
        subRowId(orig.subRowId),
        cellOffset(orig.cellOffset) {}

  RBCellCoord& operator=(const RBCellCoord& orig) {
    coreRowId = orig.coreRowId;
    subRowId = orig.subRowId;
    cellOffset = orig.cellOffset;
    return (*this);
  }

  friend std::ostream& operator<<(std::ostream& out, const RBCellCoord cc);
};

typedef RBCellCoord RBCC;

inline std::ostream& operator<<(std::ostream& out, const RBCellCoord cc) {
  out << "CellCoord {" << cc.coreRowId << "," << cc.subRowId << ","
      << cc.cellOffset << "}";
  return out;
}

#endif
