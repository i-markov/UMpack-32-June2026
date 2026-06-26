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

// Author DAP, Created Mar 10, 2005

#include "cutPosition.h"

#include <ABKCommon/abklimits.h>

#include <map>

using std::cerr;
using std::cout;
using std::endl;
using std::max;
using std::min;
using std::pair;
using std::vector;

CutLinePositionOracle::CutLinePositionOracle(const CapoBin& bin) : _bin(bin) {}

double CutLinePositionOracle::getVerticalCut(bool& certain) const {
  const RBPCoreRow& cr = _bin.getRow(0);
  double LeftMargin = 0.40;
  double RightMargin = 0.40;
  double minX = _bin.getBBox().xMin;
  double maxX = _bin.getBBox().xMax;
  if (maxX - minX <= cr.getSpacing())  // just 1 site
  {
    certain = false;
    return ((maxX + minX) / 2.0);
  }
  double minXLegal = minX + LeftMargin * _bin.getBBox().getWidth();
  // round minXLegal to next site
  double srStartMin = cr.findSubRow(minXLegal)->getXStart();
  minXLegal -= srStartMin;
  minXLegal /= cr.getSpacing();  // num sites
  minXLegal += 0.5;
  minXLegal = floor(minXLegal);
  minXLegal *= cr.getSpacing();
  minXLegal += srStartMin;
  minXLegal = max(minX + cr.getSpacing(), minXLegal);

  double maxXLegal = maxX - RightMargin * _bin.getBBox().getWidth();
  abkfatal(maxXLegal >= minXLegal, "Cutline Oracle Says: No Valid Cutlines");
  double srStartMax = cr.findSubRow(maxXLegal)->getXStart();
  maxXLegal -= srStartMax;
  maxXLegal /= cr.getSpacing();  // num sites
  maxXLegal += 0.5;
  maxXLegal = floor(maxXLegal);
  maxXLegal *= cr.getSpacing();
  maxXLegal += srStartMax;
  maxXLegal = min(maxX - cr.getSpacing(), maxXLegal);

  // Find the position of the shortest cut line length
  // In the case of a tie, align to the boundary of a feature
  // A feature is defined by a location where the cutline length changes
  double bestCutPosition = DBL_MAX;
  double bestCutLength = DBL_MAX;
  double lastCutPosition = DBL_MAX;
  double lastCutLength = DBL_MAX;
  // TODO: Stop assuming that subrows are aligned to a global grid with cell
  // width SiteSpacing Better would be to traverse site locations and jump to
  // the beginning of the next subrow when the current subrow ends
  abkfatal(_bin.rowsBegin() != _bin.rowsEnd(),
           "Cutline Oracle Says: There are no rows!");
  double SiteSpacing = cr.getSpacing();
  for (double cur_cut_position = minXLegal; cur_cut_position <= maxXLegal;
       cur_cut_position += SiteSpacing) {
    double cutLength = cutLineLengthV(cur_cut_position);
    // cerr<<"DPDEBUG: cut position "<<cur_cut_position<<" has length
    // "<<cutLength<<endl;
    if (cutLength < bestCutLength && cutLength != lastCutLength &&
        lastCutLength != DBL_MAX)  // not the very first
    {
      bestCutLength = cutLength;
      bestCutPosition = cur_cut_position;
    }

    lastCutPosition = cur_cut_position;
    lastCutLength = cutLength;
  }

  double maxCutLineLen = 0.0;
  for (itCBCoreRow rowIt = _bin.rowsBegin(); rowIt != _bin.rowsEnd(); ++rowIt)
    maxCutLineLen += (*rowIt)->getHeight();
  if (bestCutLength == maxCutLineLen || bestCutPosition == DBL_MAX) {
    // We didnt find anything special, return the middle
    bestCutPosition = (maxX - minX) / 2.0 + minX;
    // snap middle to site
    double srStartMin = cr.findSubRow(bestCutPosition)->getXStart();
    bestCutPosition -= srStartMin;
    bestCutPosition /= cr.getSpacing();  // num sites
    bestCutPosition += 0.5;
    bestCutPosition = floor(bestCutPosition);
    bestCutPosition *= cr.getSpacing();
    bestCutPosition += srStartMin;
    bestCutPosition = max(_bin.getBBox().xMin, bestCutPosition);

    certain = false;
  } else
    certain = true;

  // abkfatal(bestCutPosition != DBL_MAX, "Cutline Oracle Says: Did not find a
  // valid cut");
  if (!(bestCutPosition > minX)) {
    cerr << "Information:  bestCutPosition: " << bestCutPosition << " minX "
         << minX << endl;
    abkfatal(bestCutPosition > minX,
             "Cutline Oracle Says: Vertical cut must be *strictly* greater "
             "than minX to avoid an infinite loop");
  }

  abkfatal(bestCutPosition < maxX,
           "Cutline Oracle Says: Vertical cut must be *strictly* less than "
           "maxX to avoid an infinite loop");

  // cerr<<"DPDEBUG: Cutline Oracle Returning Vertical Cut:
  // "<<bestCutPosition<<endl; cerr<<"DPDEBUG:  which has cut length:
  // "<<bestCutLength<<" vs. max cut length: "<<maxCutLineLen<<endl;
  return bestCutPosition;
}

unsigned CutLinePositionOracle::getHorizontalCut(bool& certain) const {
  double BottomMargin = 0.40;
  double TopMargin = 0.40;
  unsigned minRowIdx = 0;
  unsigned maxRowIdx = _bin.getRows().size() - 1;
  if (maxRowIdx == 0) return 0;

  unsigned minLegalRowIdx =
      minRowIdx + static_cast<unsigned>(TopMargin * maxRowIdx);
  minLegalRowIdx = max(minLegalRowIdx, 1u);
  unsigned maxLegalRowIdx =
      maxRowIdx - static_cast<unsigned>(BottomMargin * maxRowIdx);
  abkfatal(maxLegalRowIdx >= minLegalRowIdx,
           "Cutline Oracle Says: No Valid Cutlines");

  // cerr<<"Getting horiz cut summary: "<<endl;
  // cerr<<"NumRows: "<<_bin.getRows().size()<<endl;
  // cerr<<"MinLegal: "<<minLegalRowIdx<<endl;
  // cerr<<"MaxLegal: "<<maxLegalRowIdx<<endl;

  // Find the position of the shortest cut line length
  // In the case of a tie, align to the boundary of a feature
  // A feature is defined by a location where the cutline length changes
  unsigned bestCutPosition = UINT_MAX;
  double bestCutLength = DBL_MAX;
  unsigned lastCutPosition = UINT_MAX;
  double lastCutLength = DBL_MAX;
  double maxCutLineLen = -DBL_MAX;
  abkfatal(_bin.rowsBegin() != _bin.rowsEnd(),
           "Cutline Oracle Says: There are no rows!");
  unsigned numRows = _bin.rowsEnd() - _bin.rowsBegin();
  abkfatal(numRows >= 2,
           "Cutline Oracle Says: Cannot cut a single row horizontally");
  for (unsigned cutRowIt = minLegalRowIdx; cutRowIt <= maxLegalRowIdx;
       ++cutRowIt) {
    double cutLength = cutLineLengthH(cutRowIt);
    // cerr<<"DPDEBUG: cut position "<<cur_cut_position<<" has length
    // "<<cutLength<<endl;
    if (cutLength < bestCutLength && cutLength != lastCutLength &&
        lastCutLength != DBL_MAX)  // not the very first
    {
      bestCutLength = cutLength;
      bestCutPosition = cutRowIt;
    }

    if (cutLength > maxCutLineLen) maxCutLineLen = cutLength;

    lastCutPosition = cutRowIt;
    lastCutLength = cutLength;
  }

  if (bestCutLength >= maxCutLineLen || bestCutPosition == DBL_MAX) {
    // We didnt find anything special, return the middle
    bestCutPosition = _bin.getRows().size() / 2;
    certain = false;
  } else
    certain = true;

  // abkfatal(bestCutPosition != DBL_MAX, "Cutline Oracle Says: Did not find a
  // valid cut");
  if (!(bestCutPosition > 0)) {
    // cerr<<"Information:  bestCutPosition: "<<bestCutPosition<<" minY
    // "<<minY<<endl;
    abkfatal(bestCutPosition > 0,
             "Cutline Oracle Says: Horizontal cut must be *strictly* greater "
             "than 0 to avoid an infinite loop");
  }

  if (!(bestCutPosition < numRows)) {
    // cerr<<"Information:  bestCutPosition: "<<bestCutPosition<<" minY
    // "<<minY<<endl;
    abkfatal(bestCutPosition < numRows,
             "Cutline Oracle Says: Horizontal cut must be *strictly* less than "
             "numRows to avoid an infinite loop");
  }

  // cerr<<"DPDEBUG: Cutline Oracle Returning Horizontal Cut:
  // "<<bestCutPosition<<endl; cerr<<"DPDEBUG:  which has cut length:
  // "<<bestCutLength<<" vs. max cut length: "<<maxCutLineLen<<endl;
  return bestCutPosition;
}

// Returns the number of coordinates for which there exists
// a site on both sides of the cut.  If one or both of the
// sites across the cut are missing for a coordinate, it
// does not contribute to the overall cut line length
double CutLinePositionOracle::cutLineLengthV(double cutPosition) const {
  double cutLineLength = 0.0;
  unsigned numRows = _bin.rowsEnd() - _bin.rowsBegin();
  const RBPCoreRow& cr = *(*_bin.rowsBegin());
  if (_bin.getBBox().xMax - _bin.getBBox().xMin <
      cr.getSpacing() + cr.getSiteWidth())
    // degenerate case, there is no valid cut here, length is irrelevant
    return std::numeric_limits<double>::max();

  for (unsigned row = 0; row < numRows - 1; ++row) {
    bool leftPresent = siteLeftPresent(row, cutPosition);
    bool rightPresent = siteRightPresent(row, cutPosition);
    if (leftPresent && rightPresent) {
      double rowHeight =
          _bin.getRow(row + 1).getYCoord() - _bin.getRow(row).getYCoord();
      cutLineLength += rowHeight;
    }
  }
  bool leftPresent = siteLeftPresent(numRows - 1, cutPosition);
  bool rightPresent = siteRightPresent(numRows - 1, cutPosition);
  if (leftPresent && rightPresent) {
    double rowHeight = _bin.getRow(numRows - 1).getHeight();
    cutLineLength += rowHeight;
  }
  return cutLineLength;
}

vector<double> CutLinePositionOracle::getVerticalFeatures(void) const {
  vector<pair<double, double> > featureLocations;
  const RBPCoreRow& cr = _bin.getRow(0);
  double LeftMargin = 0.40;
  double RightMargin = 0.40;
  double minX = _bin.getBBox().xMin;
  double maxX = _bin.getBBox().xMax;

  double minXLegal = minX + LeftMargin * _bin.getBBox().getWidth();
  // round minXLegal to next site
  snapToNearestSite(cr, minXLegal);
  minXLegal = max(minX + cr.getSpacing(), minXLegal);

  double maxXLegal = maxX - RightMargin * _bin.getBBox().getWidth();
  abkfatal(maxXLegal >= minXLegal, "Cutline Oracle Says: No Valid Cutlines");
  snapToNearestSite(cr, maxXLegal);
  maxXLegal = min(maxX - cr.getSpacing(), maxXLegal);

  // Find the position of the shortest cut line length
  // In the case of a tie, align to the boundary of a feature
  // A feature is defined by a location where the cutline length changes
  double lastCutPosition = DBL_MAX;
  double lastCutLength = DBL_MAX;
  // TODO: Stop assuming that subrows are aligned to a global grid with cell
  // width SiteSpacing Better would be to traverse site locations and jump to
  // the beginning of the next subrow when the current subrow ends
  abkfatal(_bin.rowsBegin() != _bin.rowsEnd(),
           "Cutline Oracle Says: There are no rows!");
  double SiteSpacing = cr.getSpacing();
  for (double cur_cut_position = minXLegal; cur_cut_position <= maxXLegal;
       cur_cut_position += SiteSpacing) {
    double cutLength = cutLineLengthV(cur_cut_position);
    if (cutLength != lastCutLength &&
        lastCutLength != DBL_MAX)  // not the very first
    {
      // featureLocations.push_back(pair<double, double>(cur_cut_position,
      // cutLength));
      featureLocations.push_back(
          pair<double, double>(cutLength, cur_cut_position));
    }
    lastCutPosition = cur_cut_position;
    lastCutLength = cutLength;
  }

  sort(featureLocations.begin(), featureLocations.end());
  vector<double> returnVal;
  for (vector<pair<double, double> >::iterator it = featureLocations.begin();
       it != featureLocations.end(); ++it)
    returnVal.push_back(it->second);
  return returnVal;
}

vector<unsigned> CutLinePositionOracle::getHorizontalFeatures(void) const {
  // vector<pair<unsigned, double> > featureLocations; // pair<cutLength,
  // location>
  vector<pair<double, unsigned> >
      featureLocations;  // pair<cutLength, location>
  double BottomMargin = 0.40;
  double TopMargin = 0.40;
  unsigned minRowIdx = 0;
  unsigned maxRowIdx = _bin.getRows().size() - 1;
  if (maxRowIdx == 0) return vector<unsigned>();

  unsigned minLegalRowIdx =
      minRowIdx + static_cast<unsigned>(TopMargin * maxRowIdx);
  minLegalRowIdx = max(minLegalRowIdx, 1u);
  unsigned maxLegalRowIdx =
      maxRowIdx - static_cast<unsigned>(BottomMargin * maxRowIdx);
  abkfatal(maxLegalRowIdx >= minLegalRowIdx,
           "Cutline Oracle Says: No Valid Cutlines");

  unsigned lastCutPosition = UINT_MAX;
  double lastCutLength = DBL_MAX;
  abkfatal(_bin.rowsBegin() != _bin.rowsEnd(),
           "Cutline Oracle Says: There are no rows!");
  unsigned numRows = _bin.rowsEnd() - _bin.rowsBegin();
  abkfatal(numRows >= 2,
           "Cutline Oracle Says: Cannot cut a single row horizontally");
  for (unsigned cutRowIt = minLegalRowIdx; cutRowIt <= maxLegalRowIdx;
       ++cutRowIt) {
    double cutLength = cutLineLengthH(cutRowIt);
    // cerr<<"DPDEBUG: cut position "<<cur_cut_position<<" has length
    // "<<cutLength<<endl;
    if (cutLength != lastCutLength &&
        lastCutLength != DBL_MAX)  // not the very first
    {
      // featureLocations.push_back( pair<unsigned, double>(cutRowIt, cutLength)
      // );
      featureLocations.push_back(pair<double, unsigned>(cutLength, cutRowIt));
    }
    lastCutPosition = cutRowIt;
    lastCutLength = cutLength;
  }

  sort(featureLocations.begin(), featureLocations.end());
  vector<unsigned> returnVal;
  for (vector<pair<double, unsigned> >::iterator it = featureLocations.begin();
       it != featureLocations.end(); ++it)
    returnVal.push_back(it->second);
  return returnVal;
}

double CutLinePositionOracle::snapToNearestSite(const RBPCoreRow& cr,
                                                double& location) const {
  double srStartMin = cr.findSubRow(location)->getXStart();
  if (location < srStartMin) return srStartMin;
  location -= srStartMin;
  location /= cr.getSpacing();  // num sites
  location += 0.5;
  location = floor(location);
  location *= cr.getSpacing();
  location += srStartMin;
  return location;
}

// Returns the number of coordinates for which there exists
// a site on both sides of the cut.  If one or both of the
// sites across the cut are missing for a coordinate, it
// does not contribute to the overall cut line length
double CutLinePositionOracle::cutLineLengthH(unsigned cutPositionRowIdx) const {
  double cutLineLength = 0.0;
  unsigned numRows = _bin.getRows().size();
  if (numRows <= 1) return std::numeric_limits<double>::max();
  const RBPCoreRow& cr = _bin.getRow(cutPositionRowIdx);
  double SiteSpacing = cr.getSpacing();
  // double SiteWidth = cr.getSiteWidth();

  double siteLoc = _bin.getBBox().xMin;
  siteLoc = snapToNearestSite(cr, siteLoc);
  siteLoc = max(_bin.getBBox().xMin, siteLoc);
  for (double siteIt = siteLoc; siteIt < _bin.getBBox().xMax;
       siteIt += SiteSpacing) {
    bool siteBelow = siteBelowPresent(cutPositionRowIdx, siteIt);
    bool siteAbove = siteAbovePresent(cutPositionRowIdx, siteIt);
    if (siteBelow && siteAbove) {
      cutLineLength += SiteSpacing;
    }
  }
  return cutLineLength;
}

// given a row index and a site coordinate,
// return whether the site left of siteCoord exists
bool CutLinePositionOracle::siteLeftPresent(unsigned rowIdx,
                                            double siteCoord) const {
  const RBPCoreRow& cr = _bin.getRow(rowIdx);
  double SiteWidth = cr.getSiteWidth();
  double SiteSpacing = cr.getSpacing();

  // make sure siteCoord is in sane range
  abkfatal(siteCoord <= _bin.getBBox().xMax,
           "Problem inputs finding a site. please report this bug");
  abkfatal(siteCoord >= _bin.getBBox().xMin,
           "Problem inputs finding a site. please report this bug");
  // abkfatal( siteCoord <= cr.getXEnd(),        "Problem inputs finding a site.
  // please report this bug" );
  if (siteCoord > cr.getXEnd() - SiteWidth + SiteSpacing) return false;
  // abkfatal( siteCoord >= cr.getXStart(),      "Problem inputs finding a site.
  // please report this bug" );
  if (siteCoord <= cr.getXStart()) return false;

  // starting with the first subrow in this bin
  // traverse subrows until siteCoord is less or equal
  // to the coordinate at end of the current subrow plus spacing.
  // or we reach the end of the bin
  itRBPSubRow srowIt = cr.subRowsBegin();
  itRBPSubRow srowEnd = cr.subRowsEnd() - 1;
  while ((srowIt != srowEnd) &&
         (siteCoord > srowIt->getXEnd() - SiteWidth + SiteSpacing))
    srowIt++;

  // Be sure that the end of the previous subrow (if it exists) is less than
  // siteCoord
  if (srowIt != cr.subRowsBegin()) {
    abkfatal(siteCoord > ((srowIt - 1)->getXEnd() - SiteWidth + SiteSpacing),
             "Problem finding a site. please report this bug");
  } else {
    abkfatal(siteCoord >= srowIt->getXStart(),
             "Problem finding a site. please report this bug");
  }

  bool insideSubrow =
      (siteCoord >= srowIt->getXStart() + SiteSpacing) &&
      (siteCoord <= srowIt->getXEnd() - SiteWidth + SiteSpacing);
  bool beforeOrBeginSubrow = siteCoord < srowIt->getXStart() + SiteSpacing;
  if (insideSubrow && !beforeOrBeginSubrow)
    return true;
  else if (beforeOrBeginSubrow && !insideSubrow)
    return false;
  else {
    cerr << "Summary of oracle:" << endl;
    cerr << "\tsiteCoord: " << siteCoord << endl;
    cerr << "\tXStart + SiteSpacing: " << (srowIt->getXStart() + SiteSpacing)
         << endl;
    cerr << "\tXEnd - SiteWidth + SiteSpacing: "
         << (srowIt->getXEnd() - SiteWidth + SiteSpacing) << endl;
    cerr << "\tXStart + SiteSpacing: " << (srowIt->getXStart() + SiteSpacing)
         << endl;
    abkfatal(0, "Problem finding a site. please report this bug");
    return true;  // to satisfy compiler
  }
}

// given a row idx and a site coordinate,
// return whether the site right of siteCoord exists
bool CutLinePositionOracle::siteRightPresent(unsigned rowIdx,
                                             double siteCoord) const {
  const RBPCoreRow& cr = _bin.getRow(rowIdx);
  double SiteWidth = cr.getSiteWidth();
  double SiteSpacing = cr.getSpacing();

  // make sure siteCoord is in sane range
  abkfatal(siteCoord <= _bin.getBBox().xMax,
           "Problem inputs finding a site, please report this bug");
  abkfatal(siteCoord >= _bin.getBBox().xMin,
           "Problem inputs finding a site, please report this bug");
  // abkfatal( siteCoord <= cr.getXEnd(),   "Problem inputs finding a site,
  // please report this bug" );
  if (siteCoord >= cr.getXEnd() - SiteWidth + SiteSpacing) return false;
  // abkfatal( siteCoord >= cr.getXStart(), "Problem inputs finding a site,
  // please report this bug" );
  if (siteCoord <= cr.getXStart() - SiteSpacing) return false;

  // starting with the last subrow in this bin
  // traverse subrows while siteCoord is less or equal
  // to the coordinate at start of the current subrow minus spacing
  // or we reach the end of the bin
  itRBPSubRow srowEnd = cr.subRowsBegin();
  itRBPSubRow srowIt = cr.subRowsEnd() - 1;
  while (srowIt != srowEnd && siteCoord <= srowIt->getXStart() - SiteSpacing)
    srowIt--;

  // Be sure that the start of the next subrow (if it exists) minus spacing is
  // greater than siteCoord
  if (srowIt != cr.subRowsEnd() - 1) {
    abkfatal(siteCoord <= (srowIt + 1)->getXStart() - SiteSpacing,
             "Problem finding a site, please report this bug");
  } else {
    abkfatal(siteCoord <= srowIt->getXEnd(),
             "Problem finding a site, please report this bug");
  }

  bool insideSubrow = (siteCoord > srowIt->getXStart() - SiteSpacing) &&
                      (siteCoord <= srowIt->getXEnd() - SiteWidth);
  bool beforeBeginSubrow = siteCoord > srowIt->getXEnd() - SiteWidth;
  if (insideSubrow && !beforeBeginSubrow)
    return true;
  else if (beforeBeginSubrow && !insideSubrow)
    return false;
  else {
    cerr << "Summary of oracle:" << endl;
    cerr << "\tsiteCoord: " << siteCoord << endl;
    cerr << "\tXStart - SiteSpacing: " << (srowIt->getXStart() - SiteSpacing)
         << endl;
    cerr << "\tXEnd - SiteWidth: " << (srowIt->getXEnd() - SiteWidth) << endl;
    cerr << "\tXEnd: " << (srowIt->getXEnd()) << endl;
    abkfatal(0, "Problem finding a site, please report this bug");
    return true;  // to guarantee a return and satisfy compiler
  }
}

// given a row idx and a site coordinate,
// return whether the site with siteCoord in the row below rowCoord exists
bool CutLinePositionOracle::siteBelowPresent(unsigned rowIdx,
                                             double siteCoord) const {
  abkfatal(rowIdx > 0, "Trying to index a row below 0, please report this bug");
  const RBPCoreRow& cr = _bin.getRow(rowIdx - 1);
  double SiteWidth = cr.getSiteWidth();
  double SiteSpacing = cr.getSpacing();

  // make sure siteCoord is in sane range
  abkfatal(siteCoord <= _bin.getBBox().xMax,
           "Problem inputs finding a site, please report this bug");
  abkfatal(siteCoord >= _bin.getBBox().xMin,
           "Problem inputs finding a site, please report this bug");
  // abkfatal( siteCoord <= cr.getXEnd(),        "Problem inputs finding a site.
  // please report this bug" );
  if (siteCoord > cr.getXEnd() - SiteWidth + SiteSpacing) return false;
  // abkfatal( siteCoord >= cr.getXStart(),      "Problem inputs finding a site.
  // please report this bug" );
  if (siteCoord < cr.getXStart()) return false;

  // if there is room for a row below rowIdx, but above
  // rowIdx-1 then there is no site below
  double thisRowCoord = _bin.getRow(rowIdx).getYCoord();
  double belowRowCoord = cr.getYCoord();
  double rowHeight = cr.getHeight();
  abkfatal(belowRowCoord <= thisRowCoord,
           "Rows are not sorted according to y-coordinate");
  if (thisRowCoord - belowRowCoord >= rowHeight * 2) return false;

  // starting with the first subrow in this row traverse subrows
  // until the end of the subrow is greater than the left edge of the bin
  // or until we reach the last subrow
  itRBPSubRow srowIt = cr.subRowsBegin();
  itRBPSubRow srowEnd = cr.subRowsEnd() - 1;
  while (srowIt != srowEnd && siteCoord > srowIt->getXEnd()) srowIt++;

  // Be sure that the end of the previous subrow (if it exists) less than
  // siteCoord
  if (srowIt != cr.subRowsBegin()) {
    abkfatal(siteCoord > (srowIt - 1)->getXEnd(),
             "Problem finding a site, please report this bug");
  } else {
    abkfatal(siteCoord <= srowIt->getXEnd(),
             "Problem finding a site, please report this bug");
  }

  bool insideSubrow =
      (siteCoord >= srowIt->getXStart()) && (siteCoord <= srowIt->getXEnd());
  bool beforeBeginSubrow = siteCoord < srowIt->getXStart();
  if (insideSubrow && !beforeBeginSubrow)
    return true;
  else if (beforeBeginSubrow && !insideSubrow)
    return false;
  else {
    cerr << "Summary of oracle:" << endl;
    cerr << "\tsiteCoord: " << siteCoord << endl;
    cerr << "\tXStart: " << (srowIt->getXStart()) << endl;
    cerr << "\tXEnd: " << (srowIt->getXEnd()) << endl;
    abkfatal(0, "Problem finding a site, please report this bug");
    return true;  // to guarantee a return and satisfy compiler
  }
}

// given a row idx and a site coordinate,
// return whether the site with siteCoord in the row above rowCoord exists
bool CutLinePositionOracle::siteAbovePresent(unsigned rowIdx,
                                             double siteCoord) const {
  const RBPCoreRow& cr = _bin.getRow(rowIdx);
  double SiteWidth = cr.getSiteWidth();
  double SiteSpacing = cr.getSpacing();

  // make sure siteCoord is in sane range
  abkfatal(siteCoord <= _bin.getBBox().xMax,
           "Problem inputs finding a site, please report this bug");
  abkfatal(siteCoord >= _bin.getBBox().xMin,
           "Problem inputs finding a site, please report this bug");
  // abkfatal( siteCoord <= cr.getXEnd(),        "Problem inputs finding a site.
  // please report this bug" );
  if (siteCoord >= cr.getXEnd() - SiteWidth + SiteSpacing) return false;
  // abkfatal( siteCoord >= cr.getXStart(),      "Problem inputs finding a site.
  // please report this bug" );
  if (siteCoord < cr.getXStart()) return false;

  // starting with the first subrow in this row traverse subrows
  // until the end of the subrow is greater than the left edge of the bin
  // or until we reach the last subrow
  itRBPSubRow srowIt = cr.subRowsBegin();
  itRBPSubRow srowEnd = cr.subRowsEnd() - 1;
  while (srowIt != srowEnd && siteCoord > srowIt->getXEnd()) srowIt++;

  // Be sure that the end of the previous subrow (if it exists) less than
  // siteCoord
  if (srowIt != cr.subRowsBegin()) {
    abkfatal(siteCoord > (srowIt - 1)->getXEnd(),
             "Problem finding a site, please report this bug");
  } else {
    abkfatal(siteCoord <= srowIt->getXEnd(),
             "Problem finding a site, please report this bug");
  }

  bool insideSubrow =
      (siteCoord >= srowIt->getXStart()) && (siteCoord <= srowIt->getXEnd());
  bool beforeBeginSubrow = siteCoord < srowIt->getXStart();
  if (insideSubrow && !beforeBeginSubrow)
    return true;
  else if (beforeBeginSubrow && !insideSubrow)
    return false;
  else {
    cerr << "Summary of oracle:" << endl;
    cerr << "\tsiteCoord: " << siteCoord << endl;
    cerr << "\tXStart: " << (srowIt->getXStart()) << endl;
    cerr << "\tXEnd: " << (srowIt->getXEnd()) << endl;
    abkfatal(0, "Problem finding a site, please report this bug");
    return true;  // to guarantee a return and satisfy compiler
  }
}
