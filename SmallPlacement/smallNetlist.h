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

// Small Placement Interface
//  - Small Netlist
//    created by Andrew Caldwell  981206

#ifndef __SMALLPL_NETLIST_H__
#define __SMALLPL_NETLIST_H__

#include <algorithm>
#include <cfloat>
#include <climits>

#include "ABKCommon/abkassert.h"
#include "Placement/placement.h"

class SmallPlacementNetlist {
  /***********************************
      Point * _pinOffsets;
          //of size #cells*#nets.  pinOffsets[cellId][netId]
          //gives the offset of the pin connecting the cell
          //and the net.  If they are not connected, the Point
          //has value (DBL_MAX, DBL_MAX).  Nets are assumed to connect
          //to cells only once. In the event of multiple net->cell
          //connections (or equiv. pins), it is suggested to use a point
          //in the middle of the pins.
   ***********************************/

  BBox* _pinOffsets;
  // aak: changed pinOffsets... it is still size of #cells * #nets.
  //      a pin's bbox should be interpreted as the box surrounding
  //      the equivalent pins; the [XMin,YMin] part of the bbox is
  //      the offset of the pin connecting the cell and the net.

  double* _minXPinOffset;
  // royj: keep the minimum pin offset per net, necessary for the
  //       proper functioning of inevitable WL calculations

  unsigned _numCells;
  unsigned _numNets;

 public:
  void calculateMinXPinOffsets(void) {
    if (_minXPinOffset == NULL) _minXPinOffset = new double[_numNets];

    for (unsigned n = 0; n < _numNets; ++n) {
      _minXPinOffset[n] = DBL_MAX;
      for (unsigned c = 0; c < _numCells; ++c) {
        if (!_pinOffsets[c * _numNets + n].isEmpty())
          _minXPinOffset[n] =
              std::min(_pinOffsets[c * _numNets + n].xMin, _minXPinOffset[n]);
      }
    }
  }

  void setPinOffset(unsigned cellId, unsigned netId, const Point& offset) {
    abkassert(_pinOffsets != 0, "must call setup before using setOffset");
    abkassert(cellId < _numCells, "cellId out of range");
    abkassert(netId < _numNets, "netId out of range");

    _pinOffsets[cellId * _numNets + netId].clear();
    _pinOffsets[cellId * _numNets + netId] += offset;

    /***********************************
            _pinOffsets[cellId*_numNets+netId].x = offset.x;
            _pinOffsets[cellId*_numNets+netId].y = offset.y;
     ***********************************/
  }

  void setPinOffset(unsigned cellId, unsigned netId, const BBox& offset) {
    abkassert(_pinOffsets != 0, "must call setup before using setOffset");
    abkassert(cellId < _numCells, "cellId out of range");
    abkassert(netId < _numNets, "netId out of range");

    _pinOffsets[cellId * _numNets + netId] = offset;
  }

  SmallPlacementNetlist(unsigned numCells, unsigned numNets)
      : _pinOffsets(NULL),
        _minXPinOffset(NULL),
        _numCells(UINT_MAX),
        _numNets(UINT_MAX) {
    setup(numCells, numNets);
  }

  SmallPlacementNetlist()
      : _pinOffsets(NULL),
        _minXPinOffset(NULL),
        _numCells(UINT_MAX),
        _numNets(UINT_MAX) {}

  void setup(unsigned numCells, unsigned numNets);

 public:
  ~SmallPlacementNetlist() {
    if (_pinOffsets != NULL) {
      delete[] _pinOffsets;
      _pinOffsets = NULL;
    }
    if (_minXPinOffset != NULL) {
      delete[] _minXPinOffset;
      _minXPinOffset = NULL;
    }
  }

  unsigned getNumCells() const { return _numCells; }
  unsigned getNumNets() const { return _numNets; }

  const BBox& getPinOffset(unsigned cellId, unsigned netId) const
  // const Point& getPinOffset(unsigned cellId, unsigned netId) const
  {
    abkassert(_pinOffsets != 0, "must call setup before using op");
    abkassert(cellId < _numCells, "cellId out of range");
    abkassert(netId < _numNets, "netId out of range");
    return _pinOffsets[cellId * _numNets + netId];
  }

  const double getMinXPinOffset(unsigned netId) const {
    abkassert(netId < _numNets, "netId out of range");
    return _minXPinOffset[netId];
  }

  void save(const char* netlistFile);

  friend std::ostream& operator<<(std::ostream& out,
                                  const SmallPlacementNetlist& netlist);
};

#endif
