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

#include <cmath>

#include "baseBinSplitter.h"
#include "capoBin.h"
#include "cutPosition.h"

using std::abs;
using std::vector;

unsigned BaseBinSplitter::shiftSplitRowAndSnap(unsigned curRow,
                                               const vector<double>& cellAreas,
                                               vector<double>& actualCaps,
                                               double localWS,
                                               double roundPct) const {
  unsigned shiftedSplitRow =
      shiftSplitRow(curRow, cellAreas, actualCaps, localWS);
  unsigned snappedSplitRow = shiftedSplitRow;
  if (_params.useCutOracle)
    snappedSplitRow = snapRowToOracleFeature(shiftedSplitRow, cellAreas,
                                             actualCaps, localWS, roundPct);
  return snappedSplitRow;
}

unsigned BaseBinSplitter::snapRowToOracleFeature(
    unsigned curRow, const vector<double>& cellAreas,
    vector<double>& actualCaps, double localWS, double roundPct) const {
  // By DAP
  // check if any of the feature locations given by the cut line
  // position oracle can satisfy whitespace constraints and
  // are located within roundPct percent away from curRow

  if ((cellAreas[0] == 0.) || (cellAreas[1] == 0.)) {
    actualCaps[0] = 0.;
    actualCaps[1] = 0.;
    for (unsigned r = 0; r < curRow; ++r)
      actualCaps[1] += _bin.getContainedAreaInCoreRow(r);
    for (unsigned r = curRow; r < _bin.getNumRows(); ++r)
      actualCaps[0] += _bin.getContainedAreaInCoreRow(r);
    return curRow;
  }

  CutLinePositionOracle cutPosition(_bin);
  vector<unsigned> features = cutPosition.getHorizontalFeatures();

  for (vector<unsigned>::iterator it = features.begin(); it != features.end();
       ++it) {
    unsigned featureRow = (*it);
    actualCaps[0] = 0.;
    actualCaps[1] = 0.;
    for (unsigned r = 0; r < featureRow; ++r)
      actualCaps[1] += _bin.getContainedAreaInCoreRow(r);
    for (unsigned r = featureRow; r < _bin.getNumRows(); ++r)
      actualCaps[0] += _bin.getContainedAreaInCoreRow(r);

    double p0WS = 1. - cellAreas[0] / actualCaps[0];
    double p1WS = 1. - cellAreas[1] / actualCaps[1];

    double curRowLoc = _bin.getRow(curRow).getYCoord();
    double featureRowLoc = _bin.getRow(featureRow).getYCoord();
    double pctDist = 100.0 * abs(curRowLoc - featureRowLoc) / _bin.getHeight();

    //        cout<<"Horizontal snapper checking row "<<featureRow<<" p0WS:
    //        "<<p0WS<<
    //            " p1WS: "<<p1WS<<" localWS: "<<localWS<<" pctDist:
    //            "<<pctDist<<endl;

    if (pctDist < roundPct && p0WS >= localWS && p1WS >= localWS) {
      //            cout<<"Horizontal snapper accepting row "<<curRow+1<<endl;
      return featureRow;
    }
  }
  //    cout<<"Horizontal snapper not making any changes"<<endl;

  actualCaps[0] = 0.;
  actualCaps[1] = 0.;
  for (unsigned r = 0; r < curRow; ++r)
    actualCaps[1] += _bin.getContainedAreaInCoreRow(r);
  for (unsigned r = curRow; r < _bin.getNumRows(); ++r)
    actualCaps[0] += _bin.getContainedAreaInCoreRow(r);
  return curRow;
}

double CapoBin::shiftXSplitAndSnap(double xsplit, double localWS, double c0,
                                   double c1, double roundPct,
                                   vector<double>& partCaps) const {
  double shiftedXSplit = shiftXSplit(xsplit, localWS, c0, c1, partCaps);
  double snappedXSplit = shiftedXSplit;
  if (_capo.getParams().splitterParams.useCutOracle)
    snappedXSplit = skylineSnapingForOracle(shiftedXSplit, localWS, c0, c1,
                                            roundPct, partCaps);
  return snappedXSplit;
}

double CapoBin::skylineSnapingForOracle(double xsplit, double localWS,
                                        double c0, double c1, double roundPct,
                                        vector<double>& partCaps) const {
  // By DAP
  // check if any of the feature locations given by the cut line
  // position oracle can satisfy whitespace constraints and
  // are located within roundPct percent away from xsplit

  double p0WS = 0.0, p1WS = 0.0;

  CutLinePositionOracle cutPosition(*this);
  vector<double> features = cutPosition.getVerticalFeatures();
  for (vector<double>::iterator it = features.begin(); it != features.end();
       ++it) {
    double featureXSplit = (*it);

    computePartAreas(featureXSplit, partCaps);
    p0WS = 1. - c0 / partCaps[0];
    p1WS = 1. - c1 / partCaps[1];

    //        cout<<"checking xSplit feature location: "<<featureXSplit<<" p0WS:
    //        "<<p0WS<<" p1WS: "<<p1WS<<" localWS: "<<localWS<<endl;
    //        cout<<"checking additional info: c0: "<<c0<<" partCaps[0]:
    //        "<<partCaps[0]<<" c1: "<<c1<<" partCaps[1]: "<<partCaps[1]<<endl;

    double distanceAsPct = abs(xsplit - featureXSplit) / getWidth() * 100.0;

    if (distanceAsPct <= roundPct && p0WS >= localWS && p1WS >= localWS) {
      //            cout<<"Accepting feature location: "<<featureXSplit<<endl;
      return featureXSplit;
    }
  }
  //    cout<<"New skyline algorithm not making changes."<<endl;
  return xsplit;
}
