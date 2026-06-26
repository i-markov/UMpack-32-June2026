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

#include "routingResources.h"

#include <algorithm>

RoutingResources::RoutingResources()
    : _numLayers(0),
      _lowestRoutingLayer(UINT_MAX),
      _isRoutingLayer(0),
      _layerPreferences(0),
      _trackPatterns(0),
      _coreRegion() {}

void RoutingResources::setNumLayers(unsigned num) {
  _numLayers = num;
  _isRoutingLayer.resize(num);
  _layerPreferences.resize(num);
}

void RoutingResources::setLowestRoutingLayer(unsigned num) {
  _lowestRoutingLayer = num;
}

void RoutingResources::setCoreRegion(const BBox &core) { _coreRegion = core; }

void RoutingResources::setIsRoutingLayer(unsigned layerNum, bool isRouting) {
  _isRoutingLayer[layerNum] = isRouting;
}

void RoutingResources::setLayerPreference(unsigned layerNum,
                                          DirPreference dir) {
  _layerPreferences[layerNum] = dir;
}

void RoutingResources::addTrackPattern(const RoutingTrackPattern &track) {
  _trackPatterns.push_back(track);
}

itTrackPattern RoutingResources::trackPatternsBegin(void) const {
  return _trackPatterns.begin();
}

itTrackPattern RoutingResources::trackPatternsEnd(void) const {
  return _trackPatterns.end();
}

std::pair<double, double> RoutingResources::getXandYRoutingTrackLengths(void) {
  std::vector<std::pair<double, unsigned> > xPosOfLongitudinalTracks,
      yPosOfLatitudinalTracks;

  for (itTrackPattern t = _trackPatterns.begin(); t != _trackPatterns.end();
       ++t) {
    if (!_isRoutingLayer[t->getLayerNum()]) continue;
    if (t->getLayerNum() == _lowestRoutingLayer) continue;
    if (t->isLongitudinal()) {
      if (_layerPreferences[t->getLayerNum()] == HORIZ) continue;
      double xcoord = t->getStart();
      for (unsigned i = 0; i < t->getNumTracks(); ++i) {
        if (_coreRegion.xMin <= xcoord && xcoord <= _coreRegion.xMax) {
          xPosOfLongitudinalTracks.push_back(
              std::make_pair(xcoord, t->getLayerNum()));
        }
        xcoord += t->getSpacing();
      }
    } else {
      if (_layerPreferences[t->getLayerNum()] == VERT) continue;
      double ycoord = t->getStart();
      for (unsigned i = 0; i < t->getNumTracks(); ++i) {
        if (_coreRegion.yMin <= ycoord && ycoord <= _coreRegion.yMax) {
          yPosOfLatitudinalTracks.push_back(
              std::make_pair(ycoord, t->getLayerNum()));
        }
        ycoord += t->getSpacing();
      }
    }
  }

  std::sort(xPosOfLongitudinalTracks.begin(), xPosOfLongitudinalTracks.end());
  std::vector<std::pair<double, unsigned> >::iterator new_end =
      unique(xPosOfLongitudinalTracks.begin(), xPosOfLongitudinalTracks.end());
  xPosOfLongitudinalTracks.erase(new_end, xPosOfLongitudinalTracks.end());
  std::sort(yPosOfLatitudinalTracks.begin(), yPosOfLatitudinalTracks.end());
  new_end =
      unique(yPosOfLatitudinalTracks.begin(), yPosOfLatitudinalTracks.end());
  yPosOfLatitudinalTracks.erase(new_end, yPosOfLatitudinalTracks.end());

  double coreWidth = _coreRegion.getWidth();
  double coreHeight = _coreRegion.getHeight();

  std::pair<double, double> len;
  len.first = coreWidth * static_cast<double>(yPosOfLatitudinalTracks.size());
  len.second =
      coreHeight * static_cast<double>(xPosOfLongitudinalTracks.size());

  return len;
}

unsigned RoutingResources::getLowestRoutingLayer(void) const {
  return _lowestRoutingLayer;
}

BBox RoutingResources::getCoreRegion(void) const { return _coreRegion; }

bool RoutingResources::isRoutingLayer(unsigned layerNum) const {
  return _isRoutingLayer[layerNum];
}

RoutingResources::DirPreference RoutingResources::getLayerPreference(
    unsigned layerNum) const {
  return _layerPreferences[layerNum];
}
