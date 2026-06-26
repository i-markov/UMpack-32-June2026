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

// Author DAP, Mar 10, 2005

#ifndef _CUT_POSITION_H_
#define _CUT_POSITION_H_

#include "capoBin.h"

class CutLinePositionOracle {
 public:
  CutLinePositionOracle(const CapoBin&);

  double getVerticalCut(bool& certain) const;
  unsigned getHorizontalCut(bool& certain) const;

  // Returns the number of coordinates for which there exists
  // a site on both sides of the cut.  If one or both of the
  // sites across the cut are missing for a coordinate, it
  // does not contribute to the overall cut line length
  double cutLineLengthV(double cutPosition) const;
  double cutLineLengthH(unsigned cutPositionRowIdx) const;

  std::vector<double> getVerticalFeatures(void) const;
  std::vector<unsigned> getHorizontalFeatures(void) const;

 private:
  // given a row coordinate and a site coordinate,
  // return whether the site in the direction specified from that
  // point exists
  bool siteLeftPresent(unsigned rowIdx, double siteCoord) const;
  bool siteRightPresent(unsigned rowIdx, double siteCoord) const;
  bool siteBelowPresent(unsigned rowIdx, double siteCoord) const;
  bool siteAbovePresent(unsigned rowIdx, double siteCoord) const;

  double snapToNearestSite(const RBPCoreRow& cr, double& location) const;

  const CapoBin& _bin;
};

#endif
