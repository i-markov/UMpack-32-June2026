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

// Created 051103 by David Papa to move all overlap removal related code and
//       definitions out of RBPlacement.h

#ifndef _OVERLAP_REM_H_
#define _OVERLAP_REM_H_

// used for removing cell overlaps
class cellMoveData {
 public:
  unsigned sourceSRIdx, sourceCell, destSRIdx, destCell;
  Point sourceNewPos, destNewPos;
  double improvement, WLcost;

  cellMoveData()
      : sourceSRIdx(UINT_MAX),
        sourceCell(UINT_MAX),
        destSRIdx(UINT_MAX),
        destCell(UINT_MAX),
        sourceNewPos(),
        destNewPos(),
        improvement(-DBL_MAX),
        WLcost(DBL_MAX){};

  cellMoveData(unsigned _sIdx, unsigned _sCell, Point &_sNewPos, unsigned _dIdx,
               unsigned _dCell, Point &_dNewPos, double _improve, double _cost)
      : sourceSRIdx(_sIdx),
        sourceCell(_sCell),
        destSRIdx(_dIdx),
        destCell(_dCell),
        sourceNewPos(_sNewPos),
        destNewPos(_dNewPos),
        improvement(_improve),
        WLcost(_cost){};

  cellMoveData(const cellMoveData &orig)
      : sourceSRIdx(orig.sourceSRIdx),
        sourceCell(orig.sourceCell),
        destSRIdx(orig.destSRIdx),
        destCell(orig.destCell),
        sourceNewPos(orig.sourceNewPos),
        destNewPos(orig.destNewPos),
        improvement(orig.improvement),
        WLcost(orig.WLcost){};
};

// used for removing macro overlaps
class macroData {
 public:
  unsigned id;
  bool fixed;
  double width, height;

  macroData(unsigned _id, bool _fixed, double _width, double _height)
      : id(_id), fixed(_fixed), width(_width), height(_height){};
  macroData(const macroData &orig)
      : id(orig.id),
        fixed(orig.fixed),
        width(orig.width),
        height(orig.height){};
};

class greedySwapData {
 public:
  unsigned c1, c2, c3;
  Point pos1, pos2, pos3;
  double improvement;

  greedySwapData()
      : c1(UINT_MAX),
        c2(UINT_MAX),
        c3(UINT_MAX),
        pos1(),
        pos2(),
        pos3(),
        improvement(-DBL_MAX){};

  greedySwapData(unsigned _c1, unsigned _c2, unsigned _c3)
      : c1(_c1),
        c2(_c2),
        c3(_c3),
        pos1(),
        pos2(),
        pos3(),
        improvement(-DBL_MAX){};

  greedySwapData(const greedySwapData &orig)
      : c1(orig.c1),
        c2(orig.c2),
        c3(orig.c3),
        pos1(orig.pos1),
        pos2(orig.pos2),
        pos3(orig.pos3),
        improvement(orig.improvement){};
};

// used to order subrows in overlap removal
class OrderBySpace {
  const std::vector<std::pair<itRBPSubRow, double> > &_haveSpace;

 public:
  OrderBySpace(const std::vector<std::pair<itRBPSubRow, double> > &haveSpace)
      : _haveSpace(haveSpace) {}
  bool operator()(unsigned i, unsigned j) {
    return _haveSpace[i].second > _haveSpace[j].second;
  }
};

// used to order macros in overlap removal
class OrderByIncreasingArea {
  const std::vector<macroData> &_macros;

 public:
  OrderByIncreasingArea(const std::vector<macroData> &macros)
      : _macros(macros) {}
  bool operator()(unsigned i, unsigned j) {
    return (_macros[i].height * _macros[i].width) <
           (_macros[j].height * _macros[j].width);
  }
};

class OrderByDecreasingArea {
  const std::vector<macroData> &_macros;

 public:
  OrderByDecreasingArea(const std::vector<macroData> &macros)
      : _macros(macros) {}
  bool operator()(unsigned i, unsigned j) {
    return (_macros[i].height * _macros[i].width) >
           (_macros[j].height * _macros[j].width);
  }
};

#endif
