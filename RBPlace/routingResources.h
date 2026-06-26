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

#ifndef _ROUTING_RESOURCES_
#define _ROUTING_RESOURCES_

#include <vector>

#include "Geoms/bbox.h"
#include "routingTrackPattern.h"

typedef std::vector<RoutingTrackPattern>::const_iterator itTrackPattern;

class RoutingResources {
 public:
  enum DirPreference { HORIZ, VERT, NOPREFERRED };

 private:
  unsigned _numLayers;
  unsigned _lowestRoutingLayer;
  std::vector<bool> _isRoutingLayer;
  std::vector<DirPreference> _layerPreferences;
  std::vector<RoutingTrackPattern> _trackPatterns;
  BBox _coreRegion;

 public:
  RoutingResources();

  void setNumLayers(unsigned num);
  void setLowestRoutingLayer(unsigned num);
  void setCoreRegion(const BBox &core);
  void setIsRoutingLayer(unsigned layerNum, bool isRouting);
  void setLayerPreference(unsigned layerNum, DirPreference dir);
  void addTrackPattern(const RoutingTrackPattern &track);
  itTrackPattern trackPatternsBegin(void) const;
  itTrackPattern trackPatternsEnd(void) const;
  unsigned getLowestRoutingLayer(void) const;
  BBox getCoreRegion(void) const;
  bool isRoutingLayer(unsigned layerNum) const;
  DirPreference getLayerPreference(unsigned layerNum) const;

  std::pair<double, double> getXandYRoutingTrackLengths(void);
};

#endif
