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

#ifndef _ROUTING_TRACK_PATTERN_
#define _ROUTING_TRACK_PATTERN_

#include <climits>

class RoutingTrackPattern {
 private:
  unsigned _layerNum;
  bool _horizontal;  // vertical if false;
                     // horizontal specifies a horizontal pattern of vertical
                     // tracks: | | | | .... like lines of longitude
  double _start;  // X if horizontal, Y if vertical
  unsigned _numTracks;
  double _spacing;

 public:
  RoutingTrackPattern()
      : _layerNum(UINT_MAX),
        _horizontal(false),
        _start(0.),
        _numTracks(0),
        _spacing(0.) {}

  RoutingTrackPattern(unsigned layer, bool horiz, double start,
                      unsigned nTracks, double spacing)
      : _layerNum(layer),
        _horizontal(horiz),
        _start(start),
        _numTracks(nTracks),
        _spacing(spacing) {}

  RoutingTrackPattern(const RoutingTrackPattern &src)
      : _layerNum(src._layerNum),
        _horizontal(src._horizontal),
        _start(src._start),
        _numTracks(src._numTracks),
        _spacing(src._spacing) {}

  unsigned getLayerNum(void) const { return _layerNum; }
  bool isHorizontal(void) const { return _horizontal; }
  bool isLongitudinal(void) const { return isHorizontal(); }
  double getStart(void) const { return _start; }
  unsigned getNumTracks(void) const { return _numTracks; }
  double getSpacing(void) const { return _spacing; }
};

#endif
