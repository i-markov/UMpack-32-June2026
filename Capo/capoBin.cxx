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

// created by Andrew Caldwell on 09/14/97  caldwell@cs.ucla.edu

// CHANGES

// aec 981124    reworked CapoBin. Removed ClNets, and changed storage of
//       cell*s to HGFixedNode Ids.  See DOCS/CapoBin.doc for
//       a description

// aec 990203    changed dbRows to RBPCoreRows

// ilm 990319    bins now maintain a reflexive 4-directed neighborhood relation
//               and print their neighbors in op<<
// ilm 020822    ported to g++ 3.1.1

#ifdef _MSC_VER
#pragma warning(disable : 4786)
#endif

#include "capoBin.h"

#include <cfloat>

#include "ABKCommon/abkassert.h"
#include "ABKCommon/abkfunc.h"
#include "Geoms/interval.h"
#include "capoPlacer.h"

using std::cout;
using std::endl;
using std::greater;
using std::less;
using std::max;
using std::min;
using std::ostream;
using std::pair;
using std::sort;
using std::string;
using std::vector;

#include <deque>
#include <functional>

CapoBin::CapoBin(const vector<unsigned>& cellIds, itCBCoreRow rowsBegin,
                 itCBCoreRow rowsEnd, const vector<CROffset>& startOffsets,
                 const vector<CROffset>& endOffsets,
                 const vector<CapoBin*>& neighbAbove,
                 const vector<CapoBin*>& neighbBelow,
                 const vector<CapoBin*>& leftNeighb,
                 const vector<CapoBin*>& rightNeighb, CapoPlacer& capo,
                 const string& name)
    : _coreRows(rowsBegin, rowsEnd),
      _startOffsets(startOffsets),
      _endOffsets(endOffsets),
      _numSites(0),
      _capacity(0.0),
      _cellIds(cellIds),
      _oracleCellIds(),
      _ecoAllowed(false),
      _totalCellArea(0.0),
      _totalMacroArea(0.0),
      _maxRecCellArea(0.0),
      _avgRowSpacing(0.0),
      _avgCellWidth(0.0),
      _avgCellHeight(0.0),
      _largestCellArea(0.0),
      _name(name),
      _capo(capo),
      _canSplitBin(false),
      _isNewBin(true), /*_savedSoln(NULL),*/
      _centerOfGravity(DBL_MAX, DBL_MAX),
      _needFP(false),
      _floorplanned(false),
      _ecoFP(false),
#ifdef USEPATOMA
      _patomaEndCase(false),
      _patomaMacrosCommitted(false),
      _patomaPacked(false),
#endif
      _neighborsAbove(neighbAbove),
      _neighborsBelow(neighbBelow),
      _leftNeighbors(leftNeighb),
      _rightNeighbors(rightNeighb) {
  setIndexToNextIndex();
  // cerr << "coreRows.size() "<<_coreRows.size()<<endl;

  abkassert(_coreRows.size() == _endOffsets.size(),
            "Mismatch CoreRows->EndOffsets");
  abkassert(_coreRows.size() == _startOffsets.size(),
            "Mismatch CoreRows->StartOffsets");
  //  abkfatal(cellIds.size() > 0, "attempt to create bin with no cells");

  _sibling = this;

  double minSiteHeight = DBL_MAX, minSiteSpacing = DBL_MAX;
  double avgSiteHeight = 0.0, avgSiteSpacing = 0.0;

  // cout << "coreRows.size() "<<_coreRows.size();
  // cout <<"  Y cood "<<_coreRows[0]->getYCoord()<<endl;

  unsigned r;
  if (_coreRows.size() > 0) {
    double Ycoo = _coreRows[0]->getYCoord();
    for (r = 1; r != _coreRows.size(); r++) {
      double newYcoo = _coreRows[r]->getYCoord();
      _avgRowSpacing += fabs(newYcoo - Ycoo);
      Ycoo = newYcoo;
    }
    if (_coreRows.size() > 1) _avgRowSpacing /= (_coreRows.size() - 1);
  }

  for (r = 0; r != _coreRows.size(); r++) {
    const RBPCoreRow& cr = *_coreRows[r];
    unsigned sites = getContainedSitesInCoreRow(r);
    _numSites += sites;
    const RBPSite& crSite = cr.getSite();
    // by sadya in SYNP. get the correct capacity of the bin
    _capacity += sites * crSite.width;
    //  _capacity += sites*              cr.getSpacing();
    //      _capacity += sites*crSite.height*cr.getSpacing();
    //  _capacity += sites*crSite.height*crSite.width;
    minSiteHeight = min(minSiteHeight, crSite.height);
    minSiteSpacing = min(minSiteSpacing, cr.getSpacing());
    avgSiteHeight += sites * crSite.height;
    avgSiteSpacing += sites * cr.getSpacing();

    double leftEdge = cr[_startOffsets[r].first].getXStart();
    leftEdge += _startOffsets[r].second * cr.getSpacing();
    double rightEdge = cr[_endOffsets[r].first].getXStart();
    rightEdge += _endOffsets[r].second * cr.getSpacing() + cr.getSiteWidth();

    _bbox += Point(leftEdge, cr.getYCoord());
    _bbox += Point(rightEdge, cr.getYCoord() + cr.getHeight());
  }

  if (_capo.getParams().cogLocation) computeCenterOfGravity();

  // abkfatal3(_numSites != 0 || _cellIds.size()==0,
  //"bin with no sites and ",_cellIds.size()," cell(s) \n");
  abkwarn(_numSites != 0 || _cellIds.size() == 0, "bin with no sites\n");
  if (_numSites != 0) {
    _capacity *= minSiteHeight;
    avgSiteHeight /= _numSites;
    avgSiteSpacing /= _numSites;
  }
  /*
     cout << " Site Height : " << avgSiteHeight << "(avg) "
          << minSiteHeight << "(min) "<< endl;
     cout << " Site Spacing: " << avgSiteSpacing<< "(avg) "
          << minSiteSpacing<< "(min) "<< endl;
  */

  const HGraphWDimensions& hgraph = capo.getNetlistHGraph();

  // compute the total, and largest, cell area
  for (unsigned c = 0; c < _cellIds.size(); c++) {
    unsigned idx = _cellIds[c];
    _avgCellWidth += hgraph.getNodeWidth(idx);
    _avgCellHeight += hgraph.getNodeHeight(idx);
    double cellArea = hgraph.getWeight(idx);
    _totalCellArea += cellArea;
    if (hgraph.isNodeMacro(idx)) _totalMacroArea += cellArea;
    _largestCellArea = max(_largestCellArea, cellArea);
  }

  if (_cellIds.size()) {
    _avgCellWidth /= _cellIds.size();
    _avgCellHeight /= _cellIds.size();
  }

  // cout << " Average cell height :" << _avgCellHeight
  //      << " width : " << _avgCellWidth << endl;

  // calculate the median cell area now
  if (_cellIds.size()) {
    vector<unsigned>::iterator medianItr =
        _cellIds.begin() + (_cellIds.size() / 2);
    nth_element(_cellIds.begin(), medianItr, _cellIds.end(),
                CompareCellIdsByArea(hgraph));
    _medianCellArea = hgraph.getWeight(*medianItr);
  } else {
    _medianCellArea = 0;
  }

  if (_capacity < _totalCellArea &&
      _capo.getParams().verb.getForMajStats() > 1) {
    cout << "Warning: bin's cell area exceeds bin's capacity by "
         << ((_totalCellArea - _capacity) / _capacity) * 100 << "%" << endl;
    cout << " Total Cell Area: " << _totalCellArea << endl;
    cout << " Capacity:        " << _capacity << endl;
  }
}

CapoBin::CapoBin(vector<CapoBin*>& srcBins, int dir)
    : _numSites(0),
      _capacity(0.0),
      _totalCellArea(0.0),
      _totalMacroArea(0.0),
      _maxRecCellArea(0.0),
      _avgRowSpacing(0.0),
      _avgCellWidth(0.0),
      _avgCellHeight(0.0),
      _largestCellArea(0),
      _name(),
      _capo(srcBins[0]->_capo),
      _canSplitBin(false),
      _isNewBin(true), /*_savedSoln(NULL),*/
      _centerOfGravity(DBL_MAX, DBL_MAX),
      _needFP(false),
      _floorplanned(false),
      _ecoFP(false)
#ifdef USEPATOMA
      ,
      _patomaEndCase(false),
      _patomaMacrosCommitted(false),
      _patomaPacked(false)
#endif
{
  setIndexToNextIndex();
  //  if(_capo.getParams().splitterParams.savePartSolnWhenUndoing)
  //     _savedSoln = new vector<unsigned>(srcBins.size(),0);

  if (dir == 0) {
    sort(srcBins.begin(), srcBins.end(), CompareSrcBinsByX());

    _ecoAllowed = true;
    for (unsigned b = 0; b < srcBins.size(); b++) {
      const CapoBin& sBin = *srcBins[b];

      if (!_bbox.isEmpty())
        abkfatal2(
            _bbox.yMin == sBin._bbox.yMin && _bbox.yMax == sBin._bbox.yMax,
            "CapoBin ctor expected row of bins \n", " with same y-projection");

      _bbox += Point(sBin._bbox.xMin, sBin._bbox.yMin);
      _bbox += Point(sBin._bbox.xMax, sBin._bbox.yMax);
      _numSites += sBin._numSites;
      _capacity += sBin._capacity;
      _totalCellArea += sBin._totalCellArea;
      _totalMacroArea += sBin._totalMacroArea;
      _largestCellArea = max(_largestCellArea, sBin._largestCellArea);
      //        if(_capo.getParams().splitterParams.savePartSolnWhenUndoing)
      //           (*_savedSoln)[b] = sBin._cellIds.size();
      _cellIds.insert(_cellIds.end(), sBin._cellIds.begin(),
                      sBin._cellIds.end());
      _oracleCellIds.insert(_oracleCellIds.end(), sBin._oracleCellIds.begin(),
                            sBin._oracleCellIds.end());
      _obstacleCellIds.insert(_obstacleCellIds.end(),
                              sBin._obstacleCellIds.begin(),
                              sBin._obstacleCellIds.end());
      _ecoAllowed = _ecoAllowed && sBin._ecoAllowed;
    }

    sort(_obstacleCellIds.begin(), _obstacleCellIds.end());
    _obstacleCellIds.erase(
        unique(_obstacleCellIds.begin(), _obstacleCellIds.end()),
        _obstacleCellIds.end());

    if (srcBins.size() == 2) {
      const CapoBin& sBin0 = *srcBins[0];
      const CapoBin& sBin1 = *srcBins[1];
      unsigned crBin0Idx = 0;
      unsigned crBin1Idx = 0;
      unsigned numBin0Rows = sBin0._coreRows.size();
      unsigned numBin1Rows = sBin1._coreRows.size();
      while (1) {
        if (crBin0Idx >= numBin0Rows && crBin1Idx >= numBin1Rows) break;
        const RBPCoreRow* cr0;
        const RBPCoreRow* cr1;
        if (crBin0Idx >= numBin0Rows) {
          cr1 = sBin1._coreRows[crBin1Idx];
          _coreRows.push_back(sBin1._coreRows[crBin1Idx]);
          _startOffsets.push_back(sBin1._startOffsets[crBin1Idx]);
          _endOffsets.push_back(sBin1._endOffsets[crBin1Idx]);
          ++crBin1Idx;
        } else if (crBin1Idx >= numBin1Rows) {
          cr0 = sBin0._coreRows[crBin0Idx];
          _coreRows.push_back(sBin0._coreRows[crBin0Idx]);
          _startOffsets.push_back(sBin0._startOffsets[crBin0Idx]);
          _endOffsets.push_back(sBin0._endOffsets[crBin0Idx]);
          ++crBin0Idx;
        } else {
          cr0 = sBin0._coreRows[crBin0Idx];
          cr1 = sBin1._coreRows[crBin1Idx];
          if (cr0->getYCoord() == cr1->getYCoord()) {
            _coreRows.push_back(sBin0._coreRows[crBin0Idx]);
            _startOffsets.push_back(sBin0._startOffsets[crBin0Idx]);
            _endOffsets.push_back(sBin1._endOffsets[crBin1Idx]);
            ++crBin0Idx;
            ++crBin1Idx;
          } else if (cr0->getYCoord() < cr1->getYCoord()) {
            _coreRows.push_back(sBin0._coreRows[crBin0Idx]);
            _startOffsets.push_back(sBin0._startOffsets[crBin0Idx]);
            _endOffsets.push_back(sBin0._endOffsets[crBin0Idx]);
            ++crBin0Idx;
          } else {
            _coreRows.push_back(sBin1._coreRows[crBin1Idx]);
            _startOffsets.push_back(sBin1._startOffsets[crBin1Idx]);
            _endOffsets.push_back(sBin1._endOffsets[crBin1Idx]);
            ++crBin1Idx;
          }
        }
      }
    } else {
      // Warning: the following code is wrong for multiple bins
      // need to implement as above. haven't done so currently for
      // lack of time
      abkfatal(
          0,
          "Internal code error #37. Please report to the UMpack maintainers.");
      _coreRows = srcBins[0]->_coreRows;
      _startOffsets = srcBins[0]->_startOffsets;
      _endOffsets = srcBins.back()->_endOffsets;
    }

    abkfatal(_coreRows.size() == _startOffsets.size() &&
                 _startOffsets.size() == _endOffsets.size(),
             "#coreRows, startOffsets and endOffsets do not match");

    //        for(unsigned b = 0; b < srcBins.size(); b++)
    //            srcBins[b]->unLinkNeighbors();

    _leftNeighbors = srcBins[0]->_leftNeighbors;
    _rightNeighbors = srcBins.back()->_rightNeighbors;
    for (unsigned b = 0; b < srcBins.size(); b++) {
      const CapoBin& sBin = *srcBins[b];
      for (unsigned i = 0; i < sBin._neighborsAbove.size(); ++i)
        _neighborsAbove.push_back(sBin._neighborsAbove[i]);
      for (unsigned i = 0; i < sBin._neighborsBelow.size(); ++i)
        _neighborsBelow.push_back(sBin._neighborsBelow[i]);
    }
    vector<CapoBin*>::iterator new_end;
    sort(_neighborsAbove.begin(), _neighborsAbove.end(), CompareBinsByIndex());
    new_end = unique(_neighborsAbove.begin(), _neighborsAbove.end());
    _neighborsAbove.erase(new_end, _neighborsAbove.end());
    sort(_neighborsBelow.begin(), _neighborsBelow.end(), CompareBinsByIndex());
    new_end = unique(_neighborsBelow.begin(), _neighborsBelow.end());
    _neighborsBelow.erase(new_end, _neighborsBelow.end());
    linkNeighbors();

    string::size_type len1 = srcBins[0]->getName().length();
    string::size_type len2 = srcBins[1]->getName().length();
    string::size_type minlen = min(len1, len2);
    string name1(srcBins[0]->getName(), 0, minlen - 3);
    string name2(srcBins[1]->getName(), 0, minlen - 3);

    if (name1 == name2) {
      _name = name1;
    } else {
      _name = "H-Row Composite Bin";
    }
  } else {
    sort(srcBins.begin(), srcBins.end(), CompareSrcBinsByY());

    _coreRows.clear();
    _startOffsets.clear();
    _endOffsets.clear();

    _ecoAllowed = true;
    for (unsigned b = 0; b < srcBins.size(); b++) {
      const CapoBin& sBin = *srcBins[b];

      if (!_bbox.isEmpty())
        abkfatal2(
            _bbox.xMin == sBin._bbox.xMin && _bbox.xMax == sBin._bbox.xMax,
            "CapoBin ctor expected column of bins \n",
            " with same x-projection");

      _bbox += Point(sBin._bbox.xMin, sBin._bbox.yMin);
      _bbox += Point(sBin._bbox.xMax, sBin._bbox.yMax);
      _numSites += sBin._numSites;
      _capacity += sBin._capacity;
      _totalCellArea += sBin._totalCellArea;
      _totalMacroArea += sBin._totalMacroArea;
      _largestCellArea = max(_largestCellArea, sBin._largestCellArea);
      //        if(_capo.getParams().splitterParams.savePartSolnWhenUndoing)
      //           (*_savedSoln)[b] = sBin._cellIds.size();
      _cellIds.insert(_cellIds.end(), sBin._cellIds.begin(),
                      sBin._cellIds.end());
      _oracleCellIds.insert(_oracleCellIds.end(), sBin._oracleCellIds.begin(),
                            sBin._oracleCellIds.end());
      _obstacleCellIds.insert(_obstacleCellIds.end(),
                              sBin._obstacleCellIds.begin(),
                              sBin._obstacleCellIds.end());
      _ecoAllowed = _ecoAllowed && sBin._ecoAllowed;

      for (unsigned i = 0; i < sBin._coreRows.size(); i++) {
        _coreRows.push_back(sBin._coreRows[i]);
        _startOffsets.push_back(sBin._startOffsets[i]);
        _endOffsets.push_back(sBin._endOffsets[i]);
      }
    }

    sort(_obstacleCellIds.begin(), _obstacleCellIds.end());
    _obstacleCellIds.erase(
        unique(_obstacleCellIds.begin(), _obstacleCellIds.end()),
        _obstacleCellIds.end());

    //    _coreRows     = srcBins[0]->_coreRows;
    //    _startOffsets = srcBins[0]->_startOffsets;
    //    _endOffsets   = srcBins[0]->_endOffsets;

    /*    abkfatal(_coreRows.size() == _startOffsets.size() &&
      _startOffsets.size() == _endOffsets.size(),
      "#coreRows, startOffsets and endOffsets do not match");*/

    //        for(unsigned b = 0; b < srcBins.size(); b++)
    //            srcBins[b]->unLinkNeighbors();

    _neighborsAbove = srcBins.back()->_neighborsAbove;
    _neighborsBelow = srcBins[0]->_neighborsBelow;
    for (unsigned b = 0; b < srcBins.size(); b++) {
      const CapoBin& sBin = *srcBins[b];
      for (unsigned i = 0; i < sBin._leftNeighbors.size(); ++i)
        _leftNeighbors.push_back(sBin._leftNeighbors[i]);
      for (unsigned i = 0; i < sBin._rightNeighbors.size(); ++i)
        _rightNeighbors.push_back(sBin._rightNeighbors[i]);
    }
    vector<CapoBin*>::iterator new_end;
    sort(_leftNeighbors.begin(), _leftNeighbors.end(), CompareBinsByIndex());
    new_end = unique(_leftNeighbors.begin(), _leftNeighbors.end());
    _leftNeighbors.erase(new_end, _leftNeighbors.end());
    sort(_rightNeighbors.begin(), _rightNeighbors.end(), CompareBinsByIndex());
    new_end = unique(_rightNeighbors.begin(), _rightNeighbors.end());
    _rightNeighbors.erase(new_end, _rightNeighbors.end());
    linkNeighbors();

    string::size_type len1 = srcBins[0]->getName().length();
    string::size_type len2 = srcBins[1]->getName().length();
    string::size_type minlen = min(len1, len2);
    string name1(srcBins[0]->getName(), 0, minlen - 3);
    string name2(srcBins[1]->getName(), 0, minlen - 3);

    if (name1 == name2) {
      _name = name1;
    } else {
      _name = "V-Row Composite Bin";
    }
  }
  _sibling = this;

  const HGraphWDimensions& hgraph = _capo.getNetlistHGraph();

  // calculate the median cell area now
  vector<unsigned>::iterator medianItr =
      _cellIds.begin() + (_cellIds.size() / 2);
  if (medianItr == _cellIds.end())
    _medianCellArea = 0;
  else {
    nth_element(_cellIds.begin(), medianItr, _cellIds.end(),
                CompareCellIdsByArea(hgraph));
    _medianCellArea = hgraph.getWeight(*medianItr);
  }

  for (unsigned c = 0; c < _cellIds.size(); c++) {
    unsigned idx = _cellIds[c];
    _avgCellWidth += hgraph.getNodeWidth(idx);
    _avgCellHeight += hgraph.getNodeHeight(idx);
  }
  if (_cellIds.size()) {
    _avgCellWidth /= _cellIds.size();
    _avgCellHeight /= _cellIds.size();
  }

  unsigned r;
  double Ycoo;
  _avgRowSpacing = 0.0;
  if (_coreRows.size() > 0) {
    Ycoo = _coreRows[0]->getYCoord();
    for (r = 1; r < _coreRows.size(); r++) {
      double newYcoo = _coreRows[r]->getYCoord();
      _avgRowSpacing += fabs(newYcoo - Ycoo);
      Ycoo = newYcoo;
    }
    if (_coreRows.size() > 1) _avgRowSpacing /= (_coreRows.size() - 1);
  }
  if (_capo.getParams().cogLocation) computeCenterOfGravity();
}

CapoBin::CapoBin(const vector<unsigned>& cellIds,
                 const vector<const RBPCoreRow*>& candidateRows,
                 const BBox& boundary, const vector<CapoBin*>& neighbAbove,
                 const vector<CapoBin*>& neighbBelow,
                 const vector<CapoBin*>& leftNeighb,
                 const vector<CapoBin*>& rightNeighb, CapoPlacer& capo,
                 const string& name)
    : _numSites(0),
      _capacity(0.0),
      _cellIds(cellIds),
      _oracleCellIds(),
      _ecoAllowed(false),
      _totalCellArea(0.0),
      _totalMacroArea(0.0),
      _maxRecCellArea(0.0),
      _avgRowSpacing(0.0),
      _avgCellWidth(0.0),
      _avgCellHeight(0.0),
      _largestCellArea(0.0),
      _name(name),
      _capo(capo),
      _canSplitBin(false),
      _isNewBin(true), /*_savedSoln(NULL),*/
      _centerOfGravity(DBL_MAX, DBL_MAX),
      _needFP(false),
      _floorplanned(false),
      _ecoFP(false),
#ifdef USEPATOMA
      _patomaEndCase(false),
      _patomaMacrosCommitted(false),
      _patomaPacked(false),
#endif
      _neighborsAbove(neighbAbove),
      _neighborsBelow(neighbBelow),
      _leftNeighbors(leftNeighb),
      _rightNeighbors(rightNeighb) {
  setIndexToNextIndex();
  _sibling = this;
  compStartAndEndOffsets(boundary, candidateRows, _coreRows, _startOffsets,
                         _endOffsets);

  /*
  for(unsigned so = 0; so < _startOffsets.size(); so++)
  {
      cout<<"("<<_startOffsets[so].first<<" - "
           <<_startOffsets[so].second<<")"<<endl;
  }
  cout<<endl<<"EndOffsets: "<<endl;
  for(so = 0; so < _endOffsets.size(); so++)
  {
      cout<<"("<<_endOffsets[so].first<<" - "
           <<_endOffsets[so].second<<")"<<endl;
  }
  */
  double minSiteHeight = DBL_MAX, minSiteSpacing = DBL_MAX;
  double avgSiteHeight = 0.0, avgSiteSpacing = 0.0;

  unsigned r;
  if (_coreRows.size() > 0) {
    double Ycoo = _coreRows[0]->getYCoord();
    for (r = 1; r != _coreRows.size(); r++) {
      double newYcoo = _coreRows[r]->getYCoord();
      _avgRowSpacing += fabs(newYcoo - Ycoo);
      Ycoo = newYcoo;
    }
    if (_coreRows.size() > 1) _avgRowSpacing /= (_coreRows.size() - 1);
  }
  for (r = 0; r != _coreRows.size(); r++) {
    const RBPCoreRow& cr = *_coreRows[r];
    unsigned sites = getContainedSitesInCoreRow(r);
    _numSites += sites;
    const RBPSite& crSite = cr.getSite();
    _capacity += sites * crSite.width;
    //  _capacity += sites*              cr.getSpacing();
    //      _capacity += sites*crSite.height*cr.getSpacing();
    //  _capacity += sites*crSite.height*crSite.width;
    minSiteHeight = min(minSiteHeight, crSite.height);
    minSiteSpacing = min(minSiteSpacing, cr.getSpacing());
    avgSiteHeight += sites * crSite.height;
    avgSiteSpacing += sites * cr.getSpacing();

    double leftEdge = cr[_startOffsets[r].first].getXStart();
    leftEdge += _startOffsets[r].second * cr.getSpacing();
    double rightEdge = cr[_endOffsets[r].first].getXStart();
    rightEdge += _endOffsets[r].second * cr.getSpacing() + cr.getSiteWidth();

    _bbox += Point(leftEdge, cr.getYCoord());
    _bbox += Point(rightEdge, cr.getYCoord() + cr.getHeight());
  }

  if (_capo.getParams().cogLocation) computeCenterOfGravity();

  if (_bbox != boundary) {
    abkwarn(0, " boundary bbox does not match computed bbox");
    abkwarn(0, " (if running from floorplanned bins, this is probably OK)");
  }

  abkfatal(_numSites != 0, "bin with no sites");

  _capacity *= minSiteHeight;
  avgSiteHeight /= _numSites;
  avgSiteSpacing /= _numSites;
  /*
     cout << " Site Height : " << avgSiteHeight << "(avg) "
          << minSiteHeight << "(min) "<< endl;
     cout << " Site Spacing: " << avgSiteSpacing<< "(avg) "
          << minSiteSpacing<< "(min) "<< endl;
  */

  const HGraphWDimensions& hgraph = _capo.getNetlistHGraph();

  // compute the total, and largest, cell area
  for (unsigned c = 0; c < _cellIds.size(); c++) {
    unsigned idx = _cellIds[c];
    _avgCellWidth += hgraph.getNodeWidth(idx);
    _avgCellHeight += hgraph.getNodeHeight(idx);
    double cellArea = hgraph.getWeight(idx);
    _totalCellArea += cellArea;
    if (hgraph.isNodeMacro(idx)) _totalMacroArea += cellArea;
    _largestCellArea = max(_largestCellArea, cellArea);
  }
  _avgCellWidth /= _cellIds.size();
  _avgCellHeight /= _cellIds.size();

  if (_capacity < _totalCellArea &&
      _capo.getParams().verb.getForMajStats() > 1) {
    cout << "Warning: bin's cell area exceeds bin's capacity by "
         << ((_totalCellArea - _capacity) / _capacity) * 100 << "%" << endl;
    cout << " Total Cell Area: " << _totalCellArea << endl;
    cout << " Capacity:        " << _capacity << endl;
  }

  // calculate the median cell area now
  vector<unsigned>::iterator medianItr =
      _cellIds.begin() + (_cellIds.size() / 2);
  nth_element(_cellIds.begin(), medianItr, _cellIds.end(),
              CompareCellIdsByArea(hgraph));
  _medianCellArea = hgraph.getWeight(*medianItr);
}

// constructor expects srcBins to be a contiguous set of
// bins, all with the same y-span.  The constructed bin
// contains all srcBins.

CapoBin::CapoBin(vector<CapoBin*>& srcBins)
    : _numSites(0),
      _capacity(0.0),
      _ecoAllowed(false),
      _totalCellArea(0.0),
      _totalMacroArea(0.0),
      _maxRecCellArea(0.0),
      _avgRowSpacing(0.0),
      _avgCellWidth(0.0),
      _avgCellHeight(0.0),
      _largestCellArea(0),
      _name(),
      _capo(srcBins[0]->_capo),
      _canSplitBin(false),
      _isNewBin(true), /*_savedSoln(NULL),*/
      _centerOfGravity(DBL_MAX, DBL_MAX),
      _needFP(false),
      _floorplanned(false),
      _ecoFP(false)
#ifdef USEPATOMA
      ,
      _patomaEndCase(false),
      _patomaMacrosCommitted(false),
      _patomaPacked(false)
#endif
{
  setIndexToNextIndex();
  sort(srcBins.begin(), srcBins.end(), CompareSrcBinsByX());

  for (unsigned b = 0; b < srcBins.size(); b++) {
    const CapoBin& sBin = *srcBins[b];

    if (!_bbox.isEmpty())
      abkfatal2(_bbox.yMin == sBin._bbox.yMin && _bbox.yMax == sBin._bbox.yMax,
                "CapoBin ctor expected row of bins \n",
                " with same y-projection");

    _bbox += Point(sBin._bbox.xMin, sBin._bbox.yMin);
    _bbox += Point(sBin._bbox.xMax, sBin._bbox.yMax);
    _numSites += sBin._numSites;
    _capacity += sBin._capacity;
    _totalCellArea += sBin._totalCellArea;
    _totalMacroArea += sBin._totalMacroArea;
    _largestCellArea = max(_largestCellArea, sBin._largestCellArea);
    _cellIds.insert(_cellIds.end(), sBin._cellIds.begin(), sBin._cellIds.end());
  }

  if (srcBins.size() == 2) {
    const CapoBin& sBin0 = *srcBins[0];
    const CapoBin& sBin1 = *srcBins[1];
    unsigned crBin0Idx = 0;
    unsigned crBin1Idx = 0;
    unsigned numBin0Rows = sBin0._coreRows.size();
    unsigned numBin1Rows = sBin1._coreRows.size();
    while (1) {
      if (crBin0Idx >= numBin0Rows && crBin1Idx >= numBin1Rows) break;
      const RBPCoreRow* cr0;
      const RBPCoreRow* cr1;
      if (crBin0Idx >= numBin0Rows) {
        cr1 = sBin1._coreRows[crBin1Idx];
        _coreRows.push_back(sBin1._coreRows[crBin1Idx]);
        _startOffsets.push_back(sBin1._startOffsets[crBin1Idx]);
        _endOffsets.push_back(sBin1._endOffsets[crBin1Idx]);
        ++crBin1Idx;
      } else if (crBin1Idx >= numBin1Rows) {
        cr0 = sBin0._coreRows[crBin0Idx];
        _coreRows.push_back(sBin0._coreRows[crBin0Idx]);
        _startOffsets.push_back(sBin0._startOffsets[crBin0Idx]);
        _endOffsets.push_back(sBin0._endOffsets[crBin0Idx]);
        ++crBin0Idx;
      } else {
        cr0 = sBin0._coreRows[crBin0Idx];
        cr1 = sBin1._coreRows[crBin1Idx];
        if (cr0->getYCoord() == cr1->getYCoord()) {
          _coreRows.push_back(sBin0._coreRows[crBin0Idx]);
          _startOffsets.push_back(sBin0._startOffsets[crBin0Idx]);
          _endOffsets.push_back(sBin1._endOffsets[crBin1Idx]);
          ++crBin0Idx;
          ++crBin1Idx;
        } else if (cr0->getYCoord() < cr1->getYCoord()) {
          _coreRows.push_back(sBin0._coreRows[crBin0Idx]);
          _startOffsets.push_back(sBin0._startOffsets[crBin0Idx]);
          _endOffsets.push_back(sBin0._endOffsets[crBin0Idx]);
          ++crBin0Idx;
        } else {
          _coreRows.push_back(sBin1._coreRows[crBin1Idx]);
          _startOffsets.push_back(sBin1._startOffsets[crBin1Idx]);
          _endOffsets.push_back(sBin1._endOffsets[crBin1Idx]);
          ++crBin1Idx;
        }
      }
    }
  } else {
    // Warning: the following code is wrong for multiple bins
    // need to implement as above. haven't done so currently for
    // lack of time
    _coreRows = srcBins[0]->_coreRows;
    _startOffsets = srcBins[0]->_startOffsets;
    _endOffsets = srcBins.back()->_endOffsets;
  }

  if (_capo.getParams().cogLocation) computeCenterOfGravity();

  abkfatal(_coreRows.size() == _startOffsets.size() &&
               _startOffsets.size() == _endOffsets.size(),
           "#coreRows, startOffsets and endOffsets do not match");

  _name = "H-Row Composite Bin";

  const HGraphWDimensions& hgraph = _capo.getNetlistHGraph();

  // calculate the median cell area now
  vector<unsigned>::iterator medianItr =
      _cellIds.begin() + (_cellIds.size() / 2);
  nth_element(_cellIds.begin(), medianItr, _cellIds.end(),
              CompareCellIdsByArea(hgraph));
  _medianCellArea = hgraph.getWeight(*medianItr);
}

// create new bin with identical information
CapoBin::CapoBin(CapoBin* src)
    : _coreRows(src->_coreRows),
      _startOffsets(src->_startOffsets),
      _endOffsets(src->_endOffsets),
      _numSites(unsigned(src->getNumSites())),
      _capacity(src->getCapacity()),
      _ecoAllowed(src->isEcoAllowed()),
      _totalCellArea(src->getTotalCellArea()),
      _totalMacroArea(src->getTotalMacroArea()),
      _maxRecCellArea(src->_maxRecCellArea),
      _avgRowSpacing(src->getAvgRowSpacing()),
      _avgCellWidth(src->getAvgCellWidth()),
      _avgCellHeight(src->getAvgCellHeight()),
      _largestCellArea(src->getLargestCellArea()),
      _medianCellArea(src->getMedianCellArea()),
      _name(src->getName()),
      _index(src->getIndex()),
      _capo(src->_capo),
      _canSplitBin(src->canSplitBin()),
      _isNewBin(src->isNewBin()), /*_savedSoln(NULL),*/
      _centerOfGravity(src->_centerOfGravity),
      _needFP(src->_needFP),
      _floorplanned(src->_floorplanned),
      _ecoFP(src->_ecoFP)
#ifdef USEPATOMA
      ,
      _patomaEndCase(src->_patomaEndCase),
      _patomaMacrosCommitted(src->_patomaMacrosCommitted),
      _patomaPacked(src->_patomaPacked)
#endif
{
  const CapoBin& sBin = *src;

  _bbox += Point(sBin._bbox.xMin, sBin._bbox.yMin);
  _bbox += Point(sBin._bbox.xMax, sBin._bbox.yMax);
  _cellIds.insert(_cellIds.end(), sBin._cellIds.begin(), sBin._cellIds.end());
  _oracleCellIds = sBin._oracleCellIds;
  _obstacleCellIds = sBin._obstacleCellIds;

  abkfatal(_coreRows.size() == _startOffsets.size() &&
               _startOffsets.size() == _endOffsets.size(),
           "#coreRows, startOffsets and endOffsets do not match");

  //   if(src->_savedSoln != NULL)
  //      _savedSoln = new vector<unsigned>(*src->_savedSoln);
}

double CapoBin::getContainedLeftEdge(unsigned r) const {
  const CROffset& startOffset = _startOffsets[r];
  const RBPCoreRow& curRow = *_coreRows[r];

  abkassert(startOffset.first < curRow.getNumSubRows(),
            "start subRow offset is too large");

  const RBPSubRow& startSR = curRow[startOffset.first];

  abkassert(startOffset.second < startSR.getNumSites(),
            "start site offset is too large");

  double leftEdge = curRow[startOffset.first].getXStart();
  leftEdge += startOffset.second * curRow.getSpacing();

  return leftEdge;
}

double CapoBin::getContainedRightEdge(unsigned r) const {
  const CROffset& endOffset = _endOffsets[r];
  const RBPCoreRow& curRow = *_coreRows[r];

  abkassert2(endOffset.first < curRow.getNumSubRows(),
             "end subRow offset is too large: ", endOffset.first);

  const RBPSubRow& endSR = curRow[endOffset.first];

  abkassert2(endOffset.second < endSR.getNumSites(),
             "end site offset is too large: ", endOffset.second);

  double rightEdge = curRow[endOffset.first].getXStart();
  rightEdge += (endOffset.second * curRow.getSpacing() + curRow.getSiteWidth());

  return rightEdge;
}

unsigned CapoBin::getContainedSitesInCoreRow(unsigned r) const {
  const CROffset& startOffset = _startOffsets[r];
  const CROffset& endOffset = _endOffsets[r];
  const RBPCoreRow& curRow = *_coreRows[r];

  abkassert2(startOffset.first < curRow.getNumSubRows(),
             "start subRow offset is too large: ", startOffset.first);
  abkassert2(endOffset.first < curRow.getNumSubRows(),
             "end subRow offset is too large: ", endOffset.first);

  const RBPSubRow& startSR = curRow[startOffset.first];
  const RBPSubRow& endSR = curRow[endOffset.first];

  abkassert2(startOffset.second < startSR.getNumSites(),
             "start site offset is too large: ", startOffset.second);
  abkassert2(endOffset.second < endSR.getNumSites(),
             "end site offset is too large: ", endOffset.second);

  unsigned numSites = 0;
  for (unsigned sr = startOffset.first; sr <= endOffset.first; sr++) {
    int sitesInThisSR = curRow[sr].getNumSites();
    if (sr == endOffset.first) sitesInThisSR = endOffset.second + 1;
    if (sr == startOffset.first) sitesInThisSR -= startOffset.second;

    numSites += sitesInThisSR;
  }
  return numSites;
}

double CapoBin::getContainedAreaInCoreRow(unsigned r) const {
  const RBPSite& crSite = _coreRows[r]->getSite();
  return getContainedSitesInCoreRow(r) * crSite.height * crSite.width;
}

double CapoBin::getContainedAreaInBBox(const BBox& rect) const {
  double containedArea = 0;
  BBox intersect = _bbox.intersect(rect);
  if (intersect.getArea() == 0) return (0);

  for (unsigned r = 0; r < getNumRows(); r++) {
    const RBPCoreRow& curRow = *_coreRows[r];
    double yCoord = curRow.getYCoord();
    double rowHeight = curRow.getHeight();
    double thisRowHeight = rowHeight;

    if ((yCoord < intersect.yMin && (yCoord + rowHeight) <= intersect.yMin) ||
        (yCoord >= intersect.yMax))
      continue;

    if ((yCoord < intersect.yMin && (yCoord + rowHeight) >= intersect.yMin))
      thisRowHeight = (yCoord + rowHeight) - intersect.yMin;
    else if (yCoord < intersect.yMax && (yCoord + rowHeight) > intersect.yMax)
      thisRowHeight = intersect.yMax - yCoord;

    const CROffset& startOffset = _startOffsets[r];
    const CROffset& endOffset = _endOffsets[r];

    abkassert2(startOffset.first < curRow.getNumSubRows(),
               "start subRow offset is too large: ", startOffset.first);
    abkassert2(endOffset.first < curRow.getNumSubRows(),
               "end subRow offset is too large: ", endOffset.first);

    const RBPSubRow& startSR = curRow[startOffset.first];
    const RBPSubRow& endSR = curRow[endOffset.first];

    abkassert2(startOffset.second < startSR.getNumSites(),
               "start site offset is too large: ", startOffset.second);
    abkassert2(endOffset.second < endSR.getNumSites(),
               "end site offset is too large: ", endOffset.second);

    // this allows for fractional sites in the region
    double numSites = 0;
    double siteSpacing = curRow.getSpacing();
    double siteWidth = curRow.getSiteWidth();

    for (unsigned sr = startOffset.first; sr <= endOffset.first; sr++) {
      double xStart = curRow[sr].getXStart();
      double xEnd = curRow[sr].getXEnd();
      if (xStart < intersect.xMin) xStart = intersect.xMin;
      if (xEnd > intersect.xMax) xEnd = intersect.xMax;

      if (xEnd <= xStart) continue;
      double sitesInThisSR = (xEnd - xStart) / siteSpacing;

      numSites += sitesInThisSR;
    }
    containedArea += numSites * thisRowHeight * siteWidth;
  }
  return (containedArea);
}

double CapoBin::findXSplit(double c0, double c1, double roundPct,
                           double p0WSWeight, double p1WSWeight,
                           vector<double>& partCaps) const {
  // it's easier to just include all begin/end of (sub) rows,
  // even if they're outside the bin...they won't have any
  // effect
  vector<double> rowEnds;
  for (unsigned c = 0; c < _coreRows.size(); ++c) {
    const RBPCoreRow& cr = *(_coreRows[c]);

    for (itRBPSubRow srIdx = cr.subRowsBegin(); srIdx != cr.subRowsEnd();
         ++srIdx) {
      if (srIdx->getXStart() > _bbox.xMin && srIdx->getXStart() < _bbox.xMax)
        rowEnds.push_back(srIdx->getXStart());
      if (srIdx->getXEnd() > _bbox.xMin && srIdx->getXEnd() < _bbox.xMax)
        rowEnds.push_back(srIdx->getXEnd());
    }
  }

  sort(rowEnds.begin(), rowEnds.end(), less<double>());
  rowEnds.erase(unique(rowEnds.begin(), rowEnds.end()), rowEnds.end());

  double xSplit =
      findXSplitUsingSkyline(c0, c1, p0WSWeight, p1WSWeight, partCaps);

  if (xSplit <= _bbox.xMin || xSplit >= _bbox.xMax) {
    double binWidth = _bbox.xMax - _bbox.xMin;
    if (xSplit <= _bbox.xMin)
      xSplit = rint((_bbox.xMin + 0.2 * binWidth));
    else
      xSplit = rint((_bbox.xMax - 0.2 * binWidth));

    if (xSplit <= _bbox.xMin || xSplit >= _bbox.xMax)
      xSplit = rint((_bbox.xMin + _bbox.xMax) / 2.0);

    computePartAreas(xSplit, partCaps);
  }
  double maxRoundDist = getWidth() * (roundPct / 100);
  double roundedSplit = DBL_MAX;
  bool foundOne = false;
  // check to see if the desired xSplit is within
  // maxRoundDist of a break. if so..take the closest one.
  // otherwise, go with the original split.
  for (unsigned s = 0; s < rowEnds.size(); s++) {
    if (rowEnds[s] >= xSplit - maxRoundDist &&
        rowEnds[s] <= xSplit + maxRoundDist) {
      roundedSplit = rowEnds[s];
      maxRoundDist = max(rowEnds[s] - xSplit, xSplit - rowEnds[s]);
      foundOne = true;
    }
  }

  if (foundOne) {
    xSplit = roundedSplit;
    // changed the splitting-line. recompute partition areas
    computePartAreas(xSplit, partCaps);
  }

  /*
   * Xiaojian: xSplit cannot be 0% or 100%
   * otherwise it causes infinite loop
   */
  if (xSplit <= _bbox.xMin || xSplit >= _bbox.xMax) {
    if (_capo.getParams().verb.getForActions() > 1) {
      cout << "xSplit is: " << xSplit << "binBBox is: " << _bbox << endl;
      abkwarn(0, "xSplit found not within bin BBox. Forcing xSplit.\n");
    }

    double binWidth = _bbox.xMax - _bbox.xMin;
    if (xSplit <= _bbox.xMin)
      xSplit = rint((_bbox.xMin + 0.3 * binWidth));
    else
      xSplit = rint((_bbox.xMax - 0.3 * binWidth));

    if (xSplit <= _bbox.xMin || xSplit >= _bbox.xMax)
      xSplit = rint((_bbox.xMin + _bbox.xMax) / 2.0);

    if (xSplit <= _bbox.xMin || xSplit >= _bbox.xMax)
      xSplit = (_bbox.xMin + _bbox.xMax) / 2.0;

    computePartAreas(xSplit, partCaps);
  }

  abkfatal(equalDouble(partCaps[0] + partCaps[1], getCapacity()),
           "partCapacities do not sum to total bin capacity");

  return xSplit;
}
// find the split which produces the best match to the given
// WS weights for each partition
// put the resulting partition site areas in partCaps
double CapoBin::findXSplitUsingSkyline(double c0, double c1, double p0WSWeight,
                                       double p1WSWeight,
                                       vector<double>& partCaps) const {
  vector<double> beginEdges;
  vector<double> endEdges;
  beginEdges.push_back(_bbox.xMax + 100);  // acts as an end point
  endEdges.push_back(_bbox.xMax + 100);    // acts as an end point
  const RBPSite& crSite = _coreRows[0]->getSite();
  double spacing = _coreRows[0]->getSpacing();

  // iterate through the contained rows, and store up the beginning
  // and end of row-intervals
  for (unsigned r = 0; r < _coreRows.size(); r++) {
    const RBPCoreRow& cr = (*_coreRows[r]);
    abkfatal(cr.getSite() == crSite, "coreRows w/ different sites");
    abkfatal(cr.getSpacing() == spacing, "coreRows w/ different spacing");

    unsigned srIdx = _startOffsets[r].first;
    for (srIdx = _startOffsets[r].first; srIdx <= _endOffsets[r].first;
         srIdx++) {
      double leftEdge = cr[srIdx].getXStart();
      if (srIdx == _startOffsets[r].first)
        leftEdge += _startOffsets[r].second * cr[srIdx].getSpacing();

      abkassert(lessOrEqualDouble(_bbox.xMin, leftEdge),
                "bbox or leftEdge is wrong");
      beginEdges.push_back(leftEdge);

      double rightEdge;
      if (srIdx == _endOffsets[r].first)
        rightEdge = cr[srIdx].getXStart() +
                    _endOffsets[r].second * cr[srIdx].getSpacing() +
                    cr[srIdx].getSiteWidth();
      else
        rightEdge = cr[srIdx].getXEnd();

      abkassert(greaterOrEqualDouble(_bbox.xMax, rightEdge),
                "bbox or rightEdge is wrong");
      endEdges.push_back(rightEdge);
    }
  }

  sort(beginEdges.begin(), beginEdges.end(), less<double>());
  sort(endEdges.begin(), endEdges.end(), less<double>());

  double siteArea = crSite.height * crSite.width;
  unsigned numRows = 0;
  unsigned beginIdx = 0;
  unsigned endIdx = 0;
  vector<double> tempVec(2, 0.);
  partCaps = tempVec;

  if (c0 == 0. || c1 == 0.) {
    double xsplit;
    if (c0 == 0.) {
      if (_capo.getParams().verb.getForMajStats() > 1 ||
          _capo.getParams().verb.getForActions() > 1)
        cout << "Part 0 of V split has no cell area, giving it one third of "
                "total area"
             << endl;
      xsplit =
          max(_bbox.xMin + spacing,
              _bbox.xMin + floor(0.3333 * (_bbox.xMax - _bbox.xMin) / spacing) *
                               spacing);
    } else {
      if (_capo.getParams().verb.getForMajStats() > 1 ||
          _capo.getParams().verb.getForActions() > 1)
        cout << "Part 1 of V split has no cell area, giving it one third of "
                "total area"
             << endl;
      xsplit =
          min(_bbox.xMax - spacing,
              _bbox.xMin + floor(0.6667 * (_bbox.xMax - _bbox.xMin) / spacing) *
                               spacing);
    }

    for (double x = _bbox.xMin; lessOrEqualDouble(x, xsplit); x += spacing) {
      while (lessThanDouble(beginEdges[beginIdx], x)) {
        numRows++;
        beginIdx++;
      }

      while (lessThanDouble(endEdges[endIdx], x)) {
        numRows--;
        endIdx++;
      }

      partCaps[0] += numRows * siteArea;
    }
    partCaps[1] = getCapacity() - partCaps[0];

    return xsplit;
  }

  double areaLeftOfX = 0;
  double bestWeightDiff = DBL_MAX;
  double bestXSplit = 0;

  for (double xSplit = _bbox.xMin; lessOrEqualDouble(xSplit, _bbox.xMax);
       xSplit += spacing) {
    while (lessThanDouble(beginEdges[beginIdx], xSplit)) {
      numRows++;
      beginIdx++;
    }

    while (lessThanDouble(endEdges[endIdx], xSplit)) {
      numRows--;
      endIdx++;
    }

    areaLeftOfX += numRows * siteArea;
    double w0 = (areaLeftOfX - c0) / areaLeftOfX;
    double w1 = (_capacity - areaLeftOfX - c1) / (_capacity - areaLeftOfX);
    double weightDiff =
        (greaterOrEqualDouble(w0, 0.) || greaterOrEqualDouble(w1, 0.))
            ? fabs(p0WSWeight * w1 - p1WSWeight * w0)
            : fabs(p0WSWeight * w0 - p1WSWeight * w1);

    if (lessThanDouble(weightDiff, bestWeightDiff)) {
      bestWeightDiff = weightDiff;
      bestXSplit = xSplit;
      partCaps[0] = areaLeftOfX;
      partCaps[1] = _capacity - areaLeftOfX;
    }
  }
  return bestXSplit;
}

double CapoBin::shiftXSplit(double xsplit, double localWS, double c0, double c1,
                            vector<double>& partCaps) const {
  double spacing = _coreRows[0]->getSpacing();
  vector<double> beginEdges;
  vector<double> endEdges;
  beginEdges.push_back(_bbox.xMax + 5 * spacing);  // acts as an end point
  endEdges.push_back(_bbox.xMax + 5 * spacing);    // acts as an end point
  const RBPSite& crSite = _coreRows[0]->getSite();

  // iterate through the contained rows, and store up the beginning
  // and end of row-intervals
  for (unsigned r = 0; r < _coreRows.size(); ++r) {
    const RBPCoreRow& cr = (*_coreRows[r]);
    abkfatal(cr.getSite() == crSite, "coreRows w/ different sites");
    abkfatal(cr.getSpacing() == spacing, "coreRows w/ different spacing");

    unsigned srIdx = _startOffsets[r].first;
    for (srIdx = _startOffsets[r].first; srIdx <= _endOffsets[r].first;
         ++srIdx) {
      double leftEdge = cr[srIdx].getXStart();
      if (srIdx == _startOffsets[r].first)
        leftEdge += _startOffsets[r].second * cr[srIdx].getSpacing();

      abkassert(lessOrEqualDouble(_bbox.xMin, leftEdge),
                "bbox or leftEdge is wrong");
      beginEdges.push_back(leftEdge);

      double rightEdge;
      if (srIdx == _endOffsets[r].first)
        rightEdge = cr[srIdx].getXStart() +
                    _endOffsets[r].second * cr[srIdx].getSpacing() +
                    cr[srIdx].getSiteWidth();
      else
        rightEdge = cr[srIdx].getXEnd();

      abkassert(greaterOrEqualDouble(_bbox.xMax, rightEdge),
                "bbox or rightEdge is wrong");
      endEdges.push_back(rightEdge);
    }
  }

  sort(beginEdges.begin(), beginEdges.end(), less<double>());
  sort(endEdges.begin(), endEdges.end(), less<double>());

  double siteArea = crSite.height * crSite.width;
  unsigned numRows = 0;
  unsigned beginIdx = 0;
  unsigned endIdx = 0;
  partCaps[0] = partCaps[1] = 0.;

  if (equalDouble(c0, 0.)) {
    if (_capo.getParams().verb.getForMajStats() > 1 ||
        _capo.getParams().verb.getForActions() > 1)
      cout << "Part 0 of V split has no cell area, giving it one third of "
              "total area"
           << endl;
    xsplit = max(_bbox.xMin + spacing,
                 _bbox.xMin + 0.3333 * (_bbox.xMax - _bbox.xMin));
  } else if (equalDouble(c1, 0.)) {
    if (_capo.getParams().verb.getForMajStats() > 1 ||
        _capo.getParams().verb.getForActions() > 1)
      cout << "Part 1 of V split has no cell area, giving it one third of "
              "total area"
           << endl;
    xsplit = min(_bbox.xMax - spacing,
                 _bbox.xMin + 0.6667 * (_bbox.xMax - _bbox.xMin));
  }

  double x;
  for (x = _bbox.xMin; lessOrEqualDouble(x, xsplit); x += spacing) {
    while (lessThanDouble(beginEdges[beginIdx], x)) {
      numRows++;
      beginIdx++;
    }

    while (lessThanDouble(endEdges[endIdx], x)) {
      numRows--;
      endIdx++;
    }

    partCaps[0] += numRows * siteArea;
  }
  partCaps[1] = getCapacity() - partCaps[0];
  xsplit = x - spacing;  // in case the given xsplit isn't site aligned

  if (equalDouble(c0, 0.) || equalDouble(c1, 0.)) {
    return xsplit;
  }

  double p0WS = 1. - (c0 / partCaps[0]);
  double p1WS = 1. - (c1 / partCaps[1]);

  if (greaterOrEqualDouble(p0WS, localWS) &&
      greaterOrEqualDouble(p1WS, localWS)) {
    return xsplit;
  }

  if (lessThanDouble(p0WS,
                     p1WS))  // increment the xsplit == add more sites to p0
  {
    if (greaterOrEqualDouble(xsplit + spacing, _bbox.xMax)) {
      if (_capo.getParams().verb.getForActions() > 1)
        abkwarn(0, "Forced to return an illegal xsplit wrt localWS");
      return xsplit;
    }
    for (xsplit = xsplit + spacing; lessOrEqualDouble(xsplit, _bbox.xMax);
         xsplit += spacing) {
      while (lessThanDouble(beginEdges[beginIdx], xsplit)) {
        numRows++;
        beginIdx++;
      }

      while (lessThanDouble(endEdges[endIdx], xsplit)) {
        numRows--;
        endIdx++;
      }

      double colArea = numRows * siteArea;
      if (lessOrEqualDouble(
              partCaps[1],
              colArea)) {  // no legal solution possible so just stop now
        if (_capo.getParams().verb.getForActions() > 1)
          abkwarn(0, "Forced to return an illegal xsplit wrt localWS");
        return (xsplit - spacing);
      }

      partCaps[0] += colArea;
      partCaps[1] -= colArea;

      p0WS = 1. - (c0 / partCaps[0]);
      p1WS = 1. - (c1 / partCaps[1]);

      if (greaterOrEqualDouble(p0WS, localWS) &&
          greaterOrEqualDouble(p1WS,
                               localWS)) {  // first legal solution, take it
        return xsplit;
      } else if (greaterOrEqualDouble(p0WS, localWS) &&
                 lessThanDouble(
                     p1WS,
                     localWS)) {  // no legal solution possible so just stop now
        if (_capo.getParams().verb.getForActions() > 1)
          abkwarn(0, "Forced to return an illegal xsplit wrt localWS");
        partCaps[0] -= colArea;
        partCaps[1] += colArea;
        return (xsplit - spacing);
      }
    }
  } else  // decrement the xsplit == add more sites to p1
  {
    if (lessOrEqualDouble(xsplit - spacing, _bbox.xMin)) {
      if (_capo.getParams().verb.getForActions() > 1)
        abkwarn(0, "Forced to return an illegal xsplit wrt localWS");
      return xsplit;
    }
    for (xsplit = xsplit - spacing; greaterOrEqualDouble(xsplit, _bbox.xMin);
         xsplit -= spacing) {
      double colArea = numRows * siteArea;

      if (lessOrEqualDouble(
              partCaps[0],
              colArea)) {  // no legal solution possible so just stop now
        if (_capo.getParams().verb.getForActions() > 1)
          abkwarn(0, "Forced to return an illegal xsplit wrt localWS");
        return (xsplit + spacing);
      }

      partCaps[0] -= colArea;
      partCaps[1] += colArea;

      while (beginIdx != 0 &&
             greaterOrEqualDouble(beginEdges[beginIdx - 1], xsplit)) {
        numRows--;
        beginIdx--;
      }

      while (endIdx != 0 &&
             greaterOrEqualDouble(endEdges[endIdx - 1], xsplit)) {
        numRows++;
        endIdx--;
      }

      p0WS = 1. - (c0 / partCaps[0]);
      p1WS = 1. - (c1 / partCaps[1]);

      if (greaterOrEqualDouble(p0WS, localWS) &&
          greaterOrEqualDouble(p1WS,
                               localWS)) {  // first legal solution, take it
        return xsplit;
      } else if (lessThanDouble(p0WS, localWS) &&
                 greaterOrEqualDouble(
                     p1WS,
                     localWS)) {  // no legal solution possible so just stop now
        if (_capo.getParams().verb.getForActions() > 1)
          abkwarn(0, "Forced to return an illegal xsplit wrt localWS");
        partCaps[0] += colArea;
        partCaps[1] -= colArea;
        return (xsplit + spacing);
      }
    }
  }
  abkfatal(0, "Couldn't return an xsplit");
  return 0.;
}

// assuming the bin is split at xSplit...put the resulting
// site areas in partCaps
void CapoBin::computePartAreas(double xSplit, vector<double>& partCaps) const {
  vector<double> beginEdges;
  vector<double> endEdges;
  beginEdges.push_back(_bbox.xMax + 100);  // acts as an end point
  endEdges.push_back(_bbox.xMax + 100);    // acts as an end point
  const RBPSite& crSite = _coreRows[0]->getSite();
  double spacing = _coreRows[0]->getSpacing();
  for (unsigned r = 0; r < _coreRows.size(); r++) {
    const RBPCoreRow& cr = (*_coreRows[r]);
    abkassert(cr.getSite() == crSite, "coreRows w/ different sites");
    abkassert(cr.getSpacing() == spacing, "coreRows w/ different spacing");

    unsigned srIdx = _startOffsets[r].first;
    for (srIdx = _startOffsets[r].first; srIdx <= _endOffsets[r].first;
         srIdx++) {
      double leftEdge = cr[srIdx].getXStart();
      if (srIdx == _startOffsets[r].first)
        leftEdge += _startOffsets[r].second * cr[srIdx].getSpacing();

      if (greaterThanDouble(_bbox.xMin, leftEdge))
        cout << _bbox.xMin << " vs " << leftEdge << endl;

      abkassert(lessOrEqualDouble(_bbox.xMin, leftEdge),
                "bbox or leftEdge is wrong");
      beginEdges.push_back(leftEdge);

      double rightEdge;
      if (srIdx == _endOffsets[r].first)
        rightEdge = cr[srIdx].getXStart() +
                    _endOffsets[r].second * cr[srIdx].getSpacing() +
                    cr[srIdx].getSiteWidth();
      else
        rightEdge = cr[srIdx].getXEnd();

      abkassert(greaterOrEqualDouble(_bbox.xMax, rightEdge),
                "bbox or rightEdge is wrong");
      endEdges.push_back(rightEdge);
    }
  }
  sort(beginEdges.begin(), beginEdges.end(), less<double>());
  sort(endEdges.begin(), endEdges.end(), less<double>());

  double siteArea = crSite.height * crSite.width;
  unsigned numRows = 0;
  double areaLeftOfX = 0;
  unsigned beginIdx = 0;
  unsigned endIdx = 0;

  for (double x = _bbox.xMin; lessOrEqualDouble(x, xSplit); x += spacing) {
    while (lessThanDouble(beginEdges[beginIdx], x)) {
      numRows++;
      beginIdx++;
    }

    while (lessThanDouble(endEdges[endIdx], x)) {
      numRows--;
      endIdx++;
    }

    areaLeftOfX += numRows * siteArea;
  }

  partCaps[0] = areaLeftOfX;
  partCaps[1] = _capacity - areaLeftOfX;
}

void CapoBin::resetCellIds(const vector<unsigned>& newCellIds) {
  abkwarn(equalDouble(_totalCellArea, _totalMacroArea) || newCellIds.size() > 0,
          "attempted reset with no cells");

  _cellIds = newCellIds;
  _totalCellArea = _totalMacroArea = _largestCellArea = 0;
  _avgCellWidth = _avgCellHeight = 0;

  const HGraphWDimensions& hgraph = _capo.getNetlistHGraph();

  for (unsigned c = 0; c < _cellIds.size(); c++) {
    abkfatal(_cellIds[c] < hgraph.getNumNodes(),
             "cellId is too large during reset");

    _totalCellArea += hgraph.getWeight(_cellIds[c]);
    if (hgraph.isNodeMacro(_cellIds[c]))
      _totalMacroArea += hgraph.getWeight(_cellIds[c]);
    _largestCellArea =
        max(_largestCellArea, double(hgraph.getWeight(_cellIds[c])));
    _avgCellWidth += hgraph.getNodeWidth(_cellIds[c]);
    _avgCellHeight += hgraph.getNodeHeight(_cellIds[c]);
  }
  if (_cellIds.size()) {
    _avgCellWidth /= _cellIds.size();
    _avgCellHeight /= _cellIds.size();
  }

  if (_capacity < _totalCellArea && _capo.getParams().verb.getForMajStats()) {
    cout << "After reset: bin's cell area exceeds bin's capacity by "
         << ((_totalCellArea - _capacity) / _capacity) * 100 << "%"
         << " after resetCellIds" << endl;
    cout << " Total Cell Area: " << _totalCellArea << endl;
    cout << " Capacity:        " << _capacity << endl;
  }
}

double CapoBin::getWidestCellWidth() const {
  const HGraphWDimensions& hgraph = _capo.getNetlistHGraph();

  double widestWidth = 0;
  for (unsigned c = 0; c < _cellIds.size(); c++)
    widestWidth =
        max(widestWidth, static_cast<double>(hgraph.getNodeWidth(_cellIds[c])));

  return widestWidth;
}

void CapoBin::linkNeighbors() {
  unsigned k;
  for (k = 0; k != _neighborsAbove.size(); k++)
    _neighborsAbove[k]->_addNeighborBelow(this);
  for (k = 0; k != _neighborsBelow.size(); k++)
    _neighborsBelow[k]->_addNeighborAbove(this);
  for (k = 0; k != _leftNeighbors.size(); k++)
    _leftNeighbors[k]->_addRightNeighbor(this);
  for (k = 0; k != _rightNeighbors.size(); k++)
    _rightNeighbors[k]->_addLeftNeighbor(this);
}

void CapoBin::unLinkNeighbors() {
  unsigned k;
  for (k = 0; k != _neighborsAbove.size(); k++)
    _neighborsAbove[k]->_removeNeighborBelow(this);
  for (k = 0; k != _neighborsBelow.size(); k++)
    _neighborsBelow[k]->_removeNeighborAbove(this);
  for (k = 0; k != _leftNeighbors.size(); k++)
    _leftNeighbors[k]->_removeRightNeighbor(this);
  for (k = 0; k != _rightNeighbors.size(); k++)
    _rightNeighbors[k]->_removeLeftNeighbor(this);
}

void CapoBin::compStartAndEndOffsets(
    const BBox& boundary, const vector<const RBPCoreRow*>& candidateRows,
    vector<const RBPCoreRow*>& coreRows, vector<CROffset>& startOffsets,
    vector<CROffset>& endOffsets) {
  for (unsigned crId = 0; crId < candidateRows.size(); crId++) {
    const RBPCoreRow& cr = *(candidateRows[crId]);
    const double spacing = cr.getSpacing();

    if (cr.getXStart() >= boundary.xMax || cr.getXEnd() <= boundary.xMin ||
        cr.getYCoord() < boundary.yMin || cr.getYCoord() >= boundary.yMax)
      continue;

    CROffset stOffset(UINT_MAX, UINT_MAX);
    CROffset endOffset(UINT_MAX, UINT_MAX);
    // find left edge
    double leftEdge = boundary.xMin;
    double rightEdge = boundary.xMax;
    itRBPSubRow srIt;

    for (srIt = cr.subRowsBegin(); srIt != cr.subRowsEnd(); srIt++) {
      if ((*srIt).getXEnd() > leftEdge &&
          (*srIt).getXStart() < rightEdge)  // 1st contained site in this sr
      {
        stOffset.first = srIt - cr.subRowsBegin();
        double delta = leftEdge - (*srIt).getXStart();

        if (delta <= 0)
          stOffset.second = 0;
        else
          stOffset.second = static_cast<unsigned>(ceil(delta / spacing));
        if (stOffset.second >= srIt->getNumSites())
          stOffset.second = srIt->getNumSites() - 1;
        break;
      }
    }
    if (stOffset.first == UINT_MAX) continue;
    // abkfatal(stOffset.first != UINT_MAX,"error: contained subrow not found -
    // StartOffsets is UINT_MAX");

    // find right edge
    for (; srIt != cr.subRowsEnd(); srIt++) {
      if ((*srIt).getXStart() >= rightEdge) break;

      endOffset.first = srIt - cr.subRowsBegin();

      double delta = (*srIt).getXEnd() - rightEdge;
      if (delta <= 0)  // contain the whole subrow
      {
        endOffset.second = (*srIt).getNumSites() - 1;
        continue;
      } else {
        unsigned excludedSites = static_cast<unsigned>(floor(delta / spacing));

        endOffset.second = ((*srIt).getNumSites() - excludedSites) - 1;
        // minus -1 because the 'end' site is to be the
        // last site contained within inside the bin
        break;
      }
    }

    if (endOffset.first == UINT_MAX) continue;
    // abkfatal(endOffset.back().first != UINT_MAX,"error: contained subrow not
    // found - EndOffsets is UINT_MAX");

    _coreRows.push_back(candidateRows[crId]);
    _startOffsets.push_back(stOffset);
    _endOffsets.push_back(endOffset);

    srIt = cr.subRowsBegin() + _startOffsets.back().first;
    if (_startOffsets.back().second >= srIt->getNumSites()) {
      cout << "start site offset is too large: " << _startOffsets.back().second
           << " vs. " << srIt->getNumSites();
      abkfatal(0, "");
    }

    srIt = cr.subRowsBegin() + _endOffsets.back().first;
    if (_endOffsets.back().second >= srIt->getNumSites()) {
      cout << "end site offset is too large: " << _endOffsets.back().second
           << " vs. " << srIt->getNumSites();
      abkfatal(0, "");
    }
  }

  /* OLD IMPLEMENTATION. TO REMOVE
    for(unsigned crId = 0; crId < candidateRows.size(); crId++)
    {
        const RBPCoreRow& cr = *(candidateRows[crId]);
        const double spacing = cr.getSpacing();

        if(cr.getXStart() >= boundary.xMax || cr.getXEnd() <= boundary.xMin ||
           cr.getYCoord() < boundary.yMin  || cr.getYCoord() >= boundary.yMax)
                continue;

        coreRows.push_back(candidateRows[crId]);
        startOffsets.push_back(CROffset(UINT_MAX,UINT_MAX));
        endOffsets.  push_back(CROffset(UINT_MAX,UINT_MAX));

                //find left edge
        double leftEdge = boundary.xMin;
        itRBPSubRow srIt;
        for(srIt = cr.subRowsBegin(); srIt != cr.subRowsEnd(); srIt++)
        {
            if((*srIt).getXEnd() >= leftEdge)  //1st contained site in this sr
            {
                startOffsets.back().first = srIt-cr.subRowsBegin();
                double delta = leftEdge-(*srIt).getXStart();

                if(delta <= 0)
                    startOffsets.back().second = 0;
                else
                    startOffsets.back().second =
                                static_cast<unsigned>(ceil(delta/spacing));
                break;
            }
        }

       abkfatal(startOffsets.back().first != UINT_MAX,
                "error: contained subrow not found - StartOffsets is UINT_MAX");

                //find right edge
        double rightEdge = boundary.xMax;
        for( ; srIt != cr.subRowsEnd(); srIt++)
        {
            if((*srIt).getXStart() >= rightEdge)
                break;

            endOffsets.back().first  = srIt-cr.subRowsBegin();

            double delta = (*srIt).getXEnd()-rightEdge;
            if(delta <= 0)      //contain the whole subrow
     {
                endOffsets.back().second = (*srIt).getNumSites()-1;
                break;
     }
            else
            {
                unsigned excludedSites =
                        static_cast<unsigned>(floor(delta/spacing));

                endOffsets.back().second =
                        ((*srIt).getNumSites()-excludedSites)-1;
                        //minus -1 because the 'end' site is to be the
                        //last site contained within inside the bin
                break;
            }
        }

        abkfatal(endOffsets.back().first != UINT_MAX,
                "error: contained subrow not found - EndOffsets is UINT_MAX");
    }
  */
}

unsigned CapoBin::getNumAdjacentCells() const {
  unsigned numAdj = 0;

  const HGraphWDimensions& hgraph = _capo.getNetlistHGraph();

  const vector<unsigned>& nodeIds = getCellIds();
  bit_vector nodes(hgraph.getNumNodes(), false);  // visited
  bit_vector nets(hgraph.getNumEdges(), false);   // insident to bin

  unsigned i;
  for (i = 0; i != nodeIds.size(); i++) nodes[nodeIds[i]] = true;
  for (i = 0; i != nodeIds.size(); i++) {
    const HGFNode& node = hgraph.getNodeByIdx(nodeIds[i]);

    for (itHGFEdgeLocal edgeIt = node.edgesBegin(); edgeIt != node.edgesEnd();
         edgeIt++) {
      const HGFEdge& edge = (**edgeIt);
      if (nets[edge.getIndex()]) continue;
      nets[edge.getIndex()] = true;
      for (itHGFNodeLocal nodeIt = edge.nodesBegin(); nodeIt != edge.nodesEnd();
           nodeIt++) {
        unsigned nodeIdx = (*nodeIt)->getIndex();
        if (nodes[nodeIdx]) continue;
        nodes[nodeIdx] = true;
        numAdj++;
      }
    }
  }
  return numAdj;
}

unsigned CapoBin::getNumAdjacentCellsWithDuplicates() const {
  unsigned numAdj = 0;

  const HGraphWDimensions& hgraph = _capo.getNetlistHGraph();

  const vector<unsigned>& nodeIds = getCellIds();
  bit_vector nodes(hgraph.getNumNodes(), false);  // visited
  bit_vector nets(hgraph.getNumEdges(), false);   // insident to bin

  unsigned i;
  for (i = 0; i != nodeIds.size(); i++) nodes[nodeIds[i]] = true;
  for (i = 0; i != nodeIds.size(); i++) {
    const HGFNode& node = hgraph.getNodeByIdx(nodeIds[i]);

    for (itHGFEdgeLocal edgeIt = node.edgesBegin(); edgeIt != node.edgesEnd();
         edgeIt++) {
      const HGFEdge& edge = (**edgeIt);
      if (nets[edge.getIndex()]) continue;
      nets[edge.getIndex()] = true;
      for (itHGFNodeLocal nodeIt = edge.nodesBegin(); nodeIt != edge.nodesEnd();
           nodeIt++) {
        unsigned nodeIdx = (*nodeIt)->getIndex();
        if (nodes[nodeIdx]) continue;
        // nodes[nodeIdx]=true;   // commented out by s.m on 01/10/01
        //  this is to allow multiple count of
        //  an adjacent cell
        numAdj++;
      }
    }
  }
  return numAdj;
}

ostream& operator<<(ostream& os, const CapoBin& cb) {
  os << "   Bin ";
  os << "\"" << cb.getName() << "\"  ";
  //  os <<  " " << setbase(16) << unsigned(&cb) << setbase(10) << "   ";
  os << "[" << cb.getIndex() << "]" << endl;

  unsigned aboveN = cb._neighborsAbove.size(),
           belowN = cb._neighborsBelow.size();
  unsigned leftN = cb._leftNeighbors.size(), rightN = cb._rightNeighbors.size();
  os << "   Neighboring binIds :  ";
  if (aboveN + belowN + leftN + rightN > 0) {
    if (aboveN + belowN + leftN + rightN > 3) os << endl;
    unsigned k;
    for (k = 0; k != aboveN; k++)
      os << "[" << cb._neighborsAbove[k]->getIndex() << "] ";
    if (aboveN) os << "(above)  ";

    for (k = 0; k != belowN; k++)
      os << "[" << cb._neighborsBelow[k]->getIndex() << "] ";
    if (belowN) os << "(below)  ";

    for (k = 0; k != leftN; k++)
      os << "[" << cb._leftNeighbors[k]->getIndex() << "] ";
    if (leftN) os << "(left)  ";

    for (k = 0; k != rightN; k++)
      os << "[" << cb._rightNeighbors[k]->getIndex() << "] ";
    if (rightN) os << "(right)  ";
  } else
    os << " none ";
  os << endl;

  os << "   BBox        : " << cb._bbox;
  os << "\n CoreRows    : " << cb.getNumRows();
  if (cb._endOffsets.size() > 0)
    os << " (First CR has "
       << 1 + (cb._endOffsets[0].first - cb._startOffsets[0].first)
       << " SubRows)";
  os << "\n Sites       : " << cb._numSites;
  os << "\n Capacity    : " << cb._capacity;
  os << "\n Cells       : " << cb._cellIds.size();
  os << "\n TotCellArea : " << cb._totalCellArea;
  os << "\n TotMacroArea: " << cb._totalMacroArea;
  os << "\n LargestCell : " << cb._largestCellArea;
  os << "\n Whitespace %: " << cb.getRelativeWhitespace() * 100.0;
  os << endl;

  return os;
}

// desiredDir = 0 == vertical Split ; desiredDir = 1 == horizontal Split
bool CapoBin::findSplitPtDir(
    bool desiredDir, vector<pair<BBox, vector<unsigned> > >& groupRegionConstr,
    bool& splitDir, double& splitPt) {
  if (getNumCells() <= 200) return (false);

  // bool needToChangeDefaultFlow = false;
  double bin_Area = _bbox.getArea();
  vector<unsigned> intersectingRegions;
  for (unsigned i = 0; i < groupRegionConstr.size(); ++i) {
    BBox& rgnBBox = groupRegionConstr[i].first;
    BBox intersection = rgnBBox.intersect(_bbox);
    double intArea = intersection.getArea();
    if ((intArea / bin_Area) > 0.9375)  // region almost covers entire bin
      return (false);
    else if ((intArea / bin_Area) > 0.0025)
      intersectingRegions.push_back(i);
    else if ((intArea / bin_Area) > 0.001)
      cout << "Ignoring Region " << i << " for alignment. Region % is "
           << (intArea / bin_Area) * 100 << endl;
  }
  if (intersectingRegions.size() == 0)  // no intersecting regions found
  {
    if (desiredDir == 0)
      splitPt = rint((_bbox.xMin + _bbox.xMax) / 2);
    else
      splitPt = rint((_bbox.yMin + _bbox.yMax) / 2);
    splitDir = desiredDir;
    cout << "No intersecting region found. Proceeding with default flow "
         << endl;
    return (false);
  }

  double binArea = _capacity;
  vector<double> potentialSplitPts;
  // desiredDir = 0 == vertical Split ; desiredDir = 1 == horizontal Split
  bool partDir = desiredDir;
  bool success = false;
  double bestSplitPt = DBL_MAX;
  double relativeDiff = DBL_MAX;
  for (unsigned num = 0; num < 6; ++num) {
    cout << "Trying to find a cutline in " << partDir << " direction " << endl;
    potentialSplitPts.clear();
    if (partDir == false) {
      for (unsigned i = 0; i < intersectingRegions.size(); ++i) {
        BBox& rgnBBox = groupRegionConstr[intersectingRegions[i]].first;
        BBox intersection = rgnBBox.intersect(_bbox);
        potentialSplitPts.push_back(intersection.xMin);
        potentialSplitPts.push_back(intersection.xMax);
      }
      bestSplitPt = DBL_MAX;
      relativeDiff = DBL_MAX;
      for (unsigned i = 0; i < potentialSplitPts.size(); ++i) {
        BBox p0BBox = _bbox;
        p0BBox.xMax = potentialSplitPts[i];
        // check if this split pt intersects any other region
        bool intersectsOtherRgns = false;
        if (num < 4) {
          for (unsigned j = 0; j < intersectingRegions.size(); ++j) {
            BBox& rgnBBox = groupRegionConstr[intersectingRegions[j]].first;
            if (potentialSplitPts[i] > rgnBBox.xMin &&
                potentialSplitPts[i] < rgnBBox.xMax) {
              double rgnSpan = rgnBBox.xMax - rgnBBox.xMin;
              double minSpan = min((potentialSplitPts[i] - rgnBBox.xMin),
                                   (rgnBBox.xMax - potentialSplitPts[i]));
              if ((num < 2 && minSpan > 0.005 * rgnSpan) ||
                  (num >= 2 && minSpan > 0.05 * rgnSpan)) {
                intersectsOtherRgns = true;
                break;
              }
            }
          }
        }
        if (!intersectsOtherRgns) {
          double p0Area = getContainedAreaInBBox(p0BBox);
          double p1Area = binArea - p0Area;
          double thisRelDiff = fabs(p0Area - p1Area) / binArea;
          double minPartArea = min(p0Area / binArea, p1Area / binArea);
          if (minPartArea > 0.01) {
            if (thisRelDiff < relativeDiff) {
              relativeDiff = thisRelDiff;
              bestSplitPt = potentialSplitPts[i];
              success = true;
            }
          }
        }
      }
    } else  // horizontal split
    {
      for (unsigned i = 0; i < intersectingRegions.size(); ++i) {
        BBox& rgnBBox = groupRegionConstr[intersectingRegions[i]].first;
        BBox intersection = rgnBBox.intersect(_bbox);
        potentialSplitPts.push_back(intersection.yMin);
        potentialSplitPts.push_back(intersection.yMax);
      }
      bestSplitPt = DBL_MAX;
      relativeDiff = DBL_MAX;
      for (unsigned i = 0; i < potentialSplitPts.size(); ++i) {
        BBox p0BBox = _bbox;
        p0BBox.yMin = potentialSplitPts[i];
        // check if this split pt intersects any other region
        bool intersectsOtherRgns = false;
        if (num < 4) {
          for (unsigned j = 0; j < intersectingRegions.size(); ++j) {
            BBox& rgnBBox = groupRegionConstr[intersectingRegions[j]].first;
            if (potentialSplitPts[i] > rgnBBox.yMin &&
                potentialSplitPts[i] < rgnBBox.yMax) {
              double rgnSpan = rgnBBox.yMax - rgnBBox.yMin;
              double minSpan = min((potentialSplitPts[i] - rgnBBox.yMin),
                                   (rgnBBox.yMax - potentialSplitPts[i]));
              if ((num < 2 && minSpan > 0.005 * rgnSpan) ||
                  (num >= 2 && minSpan > 0.05 * rgnSpan)) {
                intersectsOtherRgns = true;
                break;
              }
            }
          }
        }
        if (!intersectsOtherRgns) {
          double p0Area = getContainedAreaInBBox(p0BBox);
          double p1Area = binArea - p0Area;
          double thisRelDiff = fabs(p0Area - p1Area) / binArea;
          double minPartArea = min(p0Area / binArea, p1Area / binArea);
          if (minPartArea > 0.01) {
            if (thisRelDiff < relativeDiff) {
              relativeDiff = thisRelDiff;
              bestSplitPt = potentialSplitPts[i];
              success = true;
            }
          }
        }
      }
    }
    if (!success)
      partDir = !partDir;
    else
      break;  // done
  }

  if (success)
    cout << "Found potential cutline in " << partDir << " dir at "
         << bestSplitPt << " Checking for legality" << endl;

  if (success && partDir == 0) {
    double spacing = _coreRows[0]->getSpacing();
    bestSplitPt =
        _bbox.xMin + rint((bestSplitPt - _bbox.xMin) / spacing) * spacing;
    if (bestSplitPt <= _bbox.xMin || bestSplitPt >= _bbox.xMax) {
      success = false;
      bestSplitPt = rint((_bbox.xMin + _bbox.xMax) / 2);
    }
  } else if (success && partDir == 1) {
    // bestSplitPt = _bbox.yMin +
    // rint((bestSplitPt-_bbox.yMin)/_coreRows[0]->getHeight())*_coreRows[0]->getHeight();

    unsigned r;
    for (r = 1; r < getNumRows() - 1; r++) {
      double y0 = _coreRows[r]->getYCoord();
      double y1 = _coreRows[r + 1]->getYCoord();
      if (bestSplitPt >= y0 && bestSplitPt <= y1) {
        if ((bestSplitPt - y0) <= (y1 - bestSplitPt)) {
          bestSplitPt = y0;
        } else {
          bestSplitPt = y1;
        }
        break;
      }
    }
    if (r == (getNumRows() - 1)) {
      if (bestSplitPt >= _coreRows[r]->getYCoord() && bestSplitPt <= _bbox.yMax)
        bestSplitPt = _coreRows[r]->getYCoord();
    }
    if (bestSplitPt <= _bbox.yMin || bestSplitPt >= _bbox.yMax) {
      success = false;
      bestSplitPt = rint((_bbox.yMin + _bbox.yMax) / 2);
    }
  }

  splitPt = bestSplitPt;
  splitDir = partDir;

  if (success) {
    cout << "Split point in presence of regions is " << splitPt
         << " : Split direction is ";
    if (splitDir)
      cout << "Horizontal (1)" << endl;
    else
      cout << "Vertical (0)" << endl;
  }

  return (success);
}

#ifdef USEPATOMA
bool CapoBin::patomaCondMet() const {
  const HGraphWDimensions& hgraph = _capo.getNetlistHGraph();

  bool foundMacro = false;
  vector<unsigned>::const_iterator cellIt;
  for (vector<unsigned>::const_iterator cellIt = cellIdsBegin();
       !foundMacro && (cellIt != cellIdsEnd()); ++cellIt) {
    unsigned cellId = *cellIt;
    foundMacro = hgraph.isNodeMacro(cellId);
  }
  return !_patomaMacrosCommitted && foundMacro;
}
#endif

// desiredDir = 0 == vertical Split ; desiredDir = 1 == horizontal Split
bool CapoBin::fpCondMet(bool desiredDir, double& splitPt) {
  double maxBinDim, minBinDim;

  double macroArea = 0;
  vector<double> maxMacroAreas;
  unsigned numMacros = 0;
  unsigned numStdCells = 0;
  double maxArea1 = -1;
  double maxArea2 = -1;
  int maxAreaId1 = -1;
  int maxAreaId2 = -1;
  bool foundMacro = false;
  double binUtilization = 100 * getTotalCellArea() / getCapacity();

  const HGraphWDimensions& hgraph = _capo.getNetlistHGraph();

  for (unsigned c = 0; c < _cellIds.size(); c++) {
    unsigned idx = _cellIds[c];
    if (hgraph.isNodeMacro(idx)) {
      double cellWidth = hgraph.getNodeWidth(idx);
      double cellHeight = hgraph.getNodeHeight(idx);
      double cellArea = cellWidth * cellHeight;
      macroArea += cellArea;
      numMacros++;
      maxMacroAreas.push_back(cellArea);
      if (cellArea > maxArea1) {
        maxArea2 = maxArea1;
        maxAreaId2 = maxAreaId1;
        maxArea1 = cellArea;
        maxAreaId1 = int(idx);
      }
      foundMacro = true;
    }
  }

  numStdCells = _cellIds.size() - numMacros;

  if (!foundMacro) return false;

  if (_capo.getParams().splitterParams.eco && _ecoAllowed) {
    double macroOverlap = _capo.getRBP().calcOverlapInBBox(
        _capo.getOraclePlacement(), getBBox(), OnlyMacro);
    if (equalDouble(macroOverlap, 0.)) {
      _ecoFP = true;
      _capo.markNodesNeedFP(_cellIds);
      return true;
    }
  }

  bool forceSelectiveFP = false;
  if (numMacros >= _capo.getParams().scampiThresh) {
    // force Sel. FP when there are too many macros to floorplan
    forceSelectiveFP = true;
  }

  double childWidth = 0.0, childHeight = 0.0;
  bool needFP = false;
  if (desiredDir == 0) {
    double spacing = _coreRows[0]->getSpacing();
    splitPt = getCenter().x;
    splitPt = rint(splitPt / spacing) * spacing;
    minBinDim = min(getWidth() / 2.0, getHeight());
    maxBinDim = max(getWidth() / 2.0, getHeight());
    childWidth = getWidth() / 2.0;
    childHeight = getHeight();
  } else {
    splitPt = getCenter().y;

    /*
       unsigned r;
       for(r = 0; r < getNumRows()-1; r++)
       {
       double y0 = _coreRows[r]->getYCoord();
       double y1 = _coreRows[r+1]->getYCoord();
       if(splitPt >= y0 && splitPt <= y1)
       {
       if((splitPt-y0) <= (y1 - splitPt))
       splitPt = y0;
       else
       splitPt = y1;
       break;
       }
       }
       if(r ==  (getNumRows()-1))
       {
       if(splitPt >= _coreRows[r]->getYCoord() &&
       splitPt <= _bbox.yMax)
       splitPt = _coreRows[r]->getYCoord();
       }

       double childHeight = min(splitPt - _bbox.yMin, _bbox.yMax - splitPt);
    */

    childHeight = splitPt - _bbox.yMin;
    childWidth = getWidth();

    minBinDim = min(getWidth(), childHeight);
    maxBinDim = max(getWidth(), childHeight);
  }

  double minCurrBinDim = getWidth();
  if (minCurrBinDim > getHeight()) minCurrBinDim = getHeight();
  maxMacroAreas.clear();
  for (unsigned c = 0; c < _cellIds.size(); c++) {
    unsigned idx = _cellIds[c];
    if (hgraph.isNodeMacro(idx)) {
      double cellWidth = hgraph.getNodeWidth(idx);
      double cellHeight = hgraph.getNodeHeight(idx);
      double cellArea = cellWidth * cellHeight;
      double longerCellDim = cellWidth;
      if (longerCellDim < cellHeight) longerCellDim = cellHeight;
      if (longerCellDim > minCurrBinDim * 0.10)
        maxMacroAreas.push_back(cellArea);
    }
  }

  // maxBinDim = min(getWidth(),getHeight())/2.0;
  // maxBinDim /= 1.4;
  double tolBigMacro = maxArea1 / _totalCellArea;
  tolBigMacro = max(tolBigMacro, 0.2);
  // tolBigMacro = 0.20;
  tolBigMacro = 0.10;
  minBinDim *= (1 - tolBigMacro);
  maxBinDim *= (1 - tolBigMacro);

  bool binNeedsFP = false;

  // we now have the max dimension of the potential children bin
  for (unsigned c = 0; c < _cellIds.size(); c++) {
    unsigned idx = _cellIds[c];
    if (hgraph.isNodeMacro(idx)) {
      double cellWidth = hgraph.getNodeWidth(idx);
      double cellHeight = hgraph.getNodeHeight(idx);
      bool orient1_fits = lessOrEqualDouble(cellWidth, childWidth) &&
                          lessOrEqualDouble(cellHeight, childHeight);
      bool orient2_fits = lessOrEqualDouble(cellWidth, childHeight) &&
                          lessOrEqualDouble(cellHeight, childWidth);
      bool macro_fits = orient1_fits ||
                        (!_capo.getParams().noParquetRotation && orient2_fits);
      if (!macro_fits) {  // <aaronnn> FP if node doesn't fit into child bin
                          // given allowed rotations
        needFP = true;
        if (_capo.getParams().verb.getForActions() > 0)
          cout << endl
               << "Floorplanning " << getName() << " because of condition I "
               << endl;

        if (_capo.getParams().selectiveFloorplanning || forceSelectiveFP) {
          _capo.markNodeNeedsFP(idx);
        } else {
          binNeedsFP = true;
          _capo.markNodesNeedFP(_cellIds);
          break;
        }
      }
    }
  }

  // FP Cond II is too aggresive for prediction, it assumes that the largest
  // two macros will be in the same child bin after partitioning, this may
  // not be true.  --DAP
  bool disableFPCondII = true;
  if (!binNeedsFP &&
      !disableFPCondII)  // check if largest 2 nodes fit in the children bin
  {
    if (maxAreaId1 != -1 && maxAreaId2 != -1) {
      double cellWidth1 = hgraph.getNodeWidth(maxAreaId1);
      double cellHeight1 = hgraph.getNodeHeight(maxAreaId1);
      double cellWidth2 = hgraph.getNodeWidth(maxAreaId2);
      double cellHeight2 = hgraph.getNodeHeight(maxAreaId2);

      //    double maxMacroDim1 = max(cellWidth1,cellHeight1);
      //    double maxMacroDim2 = max(cellWidth2,cellHeight2);
      double minMacroDim1 = min(cellWidth1, cellHeight1);
      double minMacroDim2 = min(cellWidth2, cellHeight2);

      double sumOfDims;
      /*
      if(maxMacroDim1 > maxMacroDim2)
        sumOfDims = maxMacroDim1+minMacroDim2;
      else
        sumOfDims = maxMacroDim2+minMacroDim1;

      sumOfDims *= 0.95;
      */
      sumOfDims = 1.0 * (minMacroDim1 + minMacroDim2);

      if (sumOfDims > maxBinDim) {
        if (_capo.getParams().verb.getForActions() > 0)
          cout << endl
               << "Floorplanning " << getName() << " because of condition II "
               << endl;
        needFP = true;

        // <aaronnn>
        if (_capo.getParams().selectiveFloorplanning || forceSelectiveFP) {
          _capo.markNodeNeedsFP(maxAreaId1);
          _capo.markNodeNeedsFP(maxAreaId2);
        } else {
          binNeedsFP = true;
          _capo.markNodesNeedFP(_cellIds);
        }
      }
    }
  }

  if (!binNeedsFP && binUtilization > 60) {
    double macroPercent = 100 * macroArea / _totalCellArea;
    if (macroPercent > 95 && numMacros <= 20) {
      if (_capo.getParams().verb.getForActions() > 0)
        cout << endl
             << "Floorplanning " << getName() << " because of condition III "
             << endl;
      needFP = true;

      binNeedsFP = true;
      _capo.markNodesNeedFP(_cellIds);
    }
  }

  if (!binNeedsFP && binUtilization > 60) {
    double macroPercent = 100 * macroArea / _totalCellArea;
    if (macroPercent > 80 && numMacros < 10 && numStdCells < 26) {
      if (_capo.getParams().verb.getForActions() > 0)
        cout << endl
             << "Floorplanning " << getName() << " because of condition IV "
             << endl;
      needFP = true;

      binNeedsFP = true;
      _capo.markNodesNeedFP(_cellIds);
    }
  }

  if (!binNeedsFP) {
    double macroPercent = 100 * macroArea / getCapacity();
    if (macroPercent > 30 && numMacros == 1) {
      if (_capo.getParams().verb.getForActions() > 0)
        cout << endl
             << "Floorplanning " << getName() << " because of condition V "
             << endl;
      needFP = true;

      binNeedsFP = true;
      _capo.markNodesNeedFP(_cellIds);
    }
  }

  if (!binNeedsFP) {
    // There is exactly one cell and it is a macro, so floorplan
    // Post floorplanning weber place it.
    if (numMacros == 1 && _cellIds.size() == 1) {
      if (_capo.getParams().verb.getForActions() > 0) {
        cout << endl
             << "Floorplanning " << getName() << " because of condition VI "
             << endl;
      }
      needFP = true;

      binNeedsFP = true;
      _capo.markNodesNeedFP(_cellIds);
    }
  }

  sort(maxMacroAreas.begin(), maxMacroAreas.end(), greater<double>());
  if (maxMacroAreas.size() > 50) maxMacroAreas.resize(50);
  double topTenMacroArea = 0.0;
  for (unsigned i = 0; i < maxMacroAreas.size(); ++i)
    topTenMacroArea += maxMacroAreas[i];

  if (!binNeedsFP) {
    // cout<<"DPDEBUG Checking FP Cond VII"<<endl;
    // cout<<"\tArea: "<<topTenMacroArea<<endl;
    // cout<<"\tNum Macros considered: "<<maxMacroAreas.size()<<endl;
    double macroPercent = 100 * topTenMacroArea / getCapacity();
    if (macroPercent > 85 &&
        numMacros <=
            20)  // <aaronnn> FP, but only if there aren't too many nodes
    {
      if (_capo.getParams().verb.getForActions() > 0) {
        cout << endl
             << "Floorplanning " << getName() << " because of condition VII "
             << endl;
      }
      needFP = true;

      binNeedsFP = true;
      _capo.markNodesNeedFP(_cellIds);
    }
  }

  if (maxMacroAreas.size() > 20) maxMacroAreas.resize(20);
  topTenMacroArea = 0.0;
  for (unsigned i = 0; i < maxMacroAreas.size(); ++i)
    topTenMacroArea += maxMacroAreas[i];

  if (!binNeedsFP) {
    // cout<<"DPDEBUG Checking FP Cond VIII"<<endl;
    // cout<<"\tArea: "<<topTenMacroArea<<endl;
    // cout<<"\tNum Macros considered: "<<maxMacroAreas.size()<<endl;
    double macroPercent = 100 * topTenMacroArea / getCapacity();
    if (macroPercent > 60 &&
        numMacros <=
            20)  // <aaronnn> FP, but only if there aren't too many nodes
    {
      if (_capo.getParams().verb.getForActions() > 0) {
        cout << endl
             << "Floorplanning " << getName() << " because of condition VIII "
             << endl;
      }
      needFP = true;

      binNeedsFP = true;
      _capo.markNodesNeedFP(_cellIds);
    }
  }

  if (!binNeedsFP) {
    double macroPercent = 100 * macroArea / _totalCellArea;
    if (macroPercent > 95 && numMacros <= 20) {
      if (_capo.getParams().verb.getForActions() > 0) {
        cout << endl
             << "Floorplanning " << getName() << " because of condition IX "
             << endl;
      }
      needFP = true;

      binNeedsFP = true;
      _capo.markNodesNeedFP(_cellIds);
    }
  }

  return (needFP);
}

bool CapoBin::fpCondMet_Strict(void) {
  double macroArea = 0;
  unsigned numMacros = 0;
  unsigned numStdCells = 0;
  double maxArea1 = -1;
  double maxArea2 = -1;
  int maxAreaId1 = -1;
  int maxAreaId2 = -1;
  bool foundMacro = false;
  double maxBinDim = max(getWidth(), getHeight());

  const HGraphWDimensions& hgraph = _capo.getNetlistHGraph();

  for (unsigned c = 0; c < _cellIds.size(); c++) {
    unsigned idx = _cellIds[c];
    if (hgraph.isNodeMacro(idx)) {
      double cellWidth = hgraph.getNodeWidth(idx);
      double cellHeight = hgraph.getNodeHeight(idx);
      double cellArea = cellWidth * cellHeight;
      macroArea += cellArea;
      numMacros++;
      if (cellArea > maxArea1) {
        maxArea2 = maxArea1;
        maxAreaId2 = maxAreaId1;
        maxArea1 = cellArea;
        maxAreaId1 = int(idx);
      }
      foundMacro = true;
    }
  }

  numStdCells = _cellIds.size() - numMacros;

  if (!foundMacro) return false;

  bool forceSelectiveFP = false;
  if (numMacros >= _capo.getParams().scampiThresh) {
    // force Sel. FP when there are too many macros to floorplan
    forceSelectiveFP = true;
  }

  bool needFP = false;

  bool binNeedsFP = false;

  // we now have the max dimension of the potential children bin
  for (unsigned c = 0; c < _cellIds.size(); c++) {
    unsigned idx = _cellIds[c];
    if (hgraph.isNodeMacro(idx)) {
      double cellWidth = hgraph.getNodeWidth(idx);
      double cellHeight = hgraph.getNodeHeight(idx);
      // double minCellDim = min(cellWidth, cellHeight);
      // if(cellWidth > minBinDim || cellHeight > minBinDim)
      if (cellWidth > maxBinDim || cellHeight > maxBinDim)
      // if(minCellDim > minBinDim)
      {
        needFP = true;
        if (_capo.getParams().verb.getForActions() > 0)
          cout << endl
               << "Backtracking from newly created " << getName()
               << " because of condition I " << endl;

        if (_capo.getParams().selectiveFloorplanning || forceSelectiveFP) {
          _capo.markNodeNeedsFP(idx);
        } else {
          binNeedsFP = true;
          _capo.markNodesNeedFP(_cellIds);
          break;
        }
        // break; <aaronnn> don't stop, find all nodes that need FPing
      }
    }
  }

  // if(!needFP) //check if largest 2 bins fit in the children bin
  if (!binNeedsFP) {
    if (maxAreaId1 != -1 && maxAreaId2 != -1) {
      double cellWidth1 = hgraph.getNodeWidth(maxAreaId1);
      double cellHeight1 = hgraph.getNodeHeight(maxAreaId1);
      double cellWidth2 = hgraph.getNodeWidth(maxAreaId2);
      double cellHeight2 = hgraph.getNodeHeight(maxAreaId2);

      //    double maxMacroDim1 = max(cellWidth1,cellHeight1);
      //    double maxMacroDim2 = max(cellWidth2,cellHeight2);
      double minMacroDim1 = min(cellWidth1, cellHeight1);
      double minMacroDim2 = min(cellWidth2, cellHeight2);

      double sumOfDims;
      /*
      if(maxMacroDim1 > maxMacroDim2)
        sumOfDims = maxMacroDim1+minMacroDim2;
      else
        sumOfDims = maxMacroDim2+minMacroDim1;

      sumOfDims *= 0.95;
      */
      sumOfDims = 1.0 * (minMacroDim1 + minMacroDim2);

      if (sumOfDims > maxBinDim) {
        if (_capo.getParams().verb.getForActions() > 0)
          cout << endl
               << "Backtracking from newly created " << getName()
               << " because of condition II " << endl;
        needFP = true;
        if (_capo.getParams().selectiveFloorplanning || forceSelectiveFP) {
          _capo.markNodeNeedsFP(maxAreaId1);
          _capo.markNodeNeedsFP(maxAreaId2);
        } else {
          binNeedsFP = true;
          _capo.markNodesNeedFP(_cellIds);
        }
      }
    }
  }

  return (needFP);
}

double CapoBin::largestMacroArea() {
  const HGraphWDimensions& hgraph = _capo.getNetlistHGraph();

  double maxArea = 0;
  for (unsigned c = 0; c < _cellIds.size(); c++) {
    unsigned idx = _cellIds[c];
    if (hgraph.isNodeMacro(idx)) {
      double cellWidth = hgraph.getNodeWidth(idx);
      double cellHeight = hgraph.getNodeHeight(idx);
      double cellArea = cellWidth * cellHeight;
      if (cellArea > maxArea) maxArea = cellArea;
    }
  }
  return (maxArea);
}

void CapoBin::removeMacros() {
  const HGraphWDimensions& hgraph = _capo.getNetlistHGraph();

  vector<unsigned> newCellIds;
  newCellIds.reserve(_cellIds.size());
  _obstacleCellIds.reserve(_obstacleCellIds.size() + _cellIds.size());
  for (unsigned i = 0; i < _cellIds.size(); ++i) {
    // <aaronnn> macros removed from bin become obstacles
    _obstacleCellIds.push_back(_cellIds[i]);

    if (!hgraph.isNodeMacro(_cellIds[i])) newCellIds.push_back(_cellIds[i]);
  }
  resetCellIds(newCellIds);
}

void CapoBin::removeFloorplannedMacros() {
  const HGraphWDimensions& hgraph = _capo.getNetlistHGraph();

  vector<unsigned> newCellIds;
  newCellIds.reserve(_cellIds.size());
  _obstacleCellIds.reserve(_obstacleCellIds.size() + _cellIds.size());
  for (unsigned i = 0; i < _cellIds.size(); ++i) {
    if (hgraph.isNodeMacro(_cellIds[i]) && _capo.nodeNeedsFP(_cellIds[i])) {
      // - macros singled-out as needing FPing,
      // - macros that are in bins that needed FPing

      // the sites for these macros would already have been removed
      // -- they become obstacles
      _obstacleCellIds.push_back(_cellIds[i]);
    } else {
      // these are the cells we want to keep in the bin
      newCellIds.push_back(_cellIds[i]);
    }
  }
  // if(newCellIds.size() > 0)
  resetCellIds(newCellIds);
}

void CapoBin::resetRows() {
  vector<const RBPCoreRow*> candidateRows;

  double rowHeight = _capo.getRBP().getNumCoreRows() > 0
                         ? _capo.getRBP().getCoreRow(0).getHeight()
                         : 1.;

  Point botPt = _bbox.getBottomLeft();
  Point topPt = _bbox.getTopRight();

  unsigned startIdx = _capo.getRBP().findCoreRowIdx(botPt);
  while (startIdx == UINT_MAX && botPt.y <= topPt.y) {
    botPt.y += rowHeight;
    startIdx = _capo.getRBP().findCoreRowIdx(botPt);
  }

  unsigned endIdx = _capo.getRBP().findCoreRowIdx(topPt);
  while (endIdx == UINT_MAX && botPt.y <= topPt.y) {
    topPt.y -= rowHeight;
    endIdx = _capo.getRBP().findCoreRowIdx(topPt);
  }

  for (unsigned i = startIdx; i <= endIdx; ++i) {
    const RBPCoreRow* cand = &_capo.getRBP().getCoreRow(i);
    candidateRows.push_back(cand);
  }

  _coreRows.clear();
  _startOffsets.clear();
  _endOffsets.clear();

  compStartAndEndOffsets(_bbox, candidateRows, _coreRows, _startOffsets,
                         _endOffsets);

  double minSiteHeight = DBL_MAX, minSiteSpacing = DBL_MAX;
  double avgSiteHeight = 0.0, avgSiteSpacing = 0.0;

  unsigned r;
  double Ycoo;
  _avgRowSpacing = 0.0;
  if (_coreRows.size() > 0) {
    Ycoo = _coreRows[0]->getYCoord();
    for (r = 1; r < _coreRows.size(); r++) {
      double newYcoo = _coreRows[r]->getYCoord();
      _avgRowSpacing += fabs(newYcoo - Ycoo);
      Ycoo = newYcoo;
    }
    if (_coreRows.size() > 1) _avgRowSpacing /= (_coreRows.size() - 1);
  }

  _numSites = 0;
  _capacity = 0;

  if (_coreRows.size() > 0) {
    _bbox.clear();

    for (r = 0; r < _coreRows.size(); r++) {
      const RBPCoreRow& cr = *_coreRows[r];
      unsigned sites = getContainedSitesInCoreRow(r);
      _numSites += sites;
      const RBPSite& crSite = cr.getSite();
      _capacity += sites * crSite.width;
      //      _capacity += sites*              cr.getSpacing();
      //      _capacity += sites*crSite.height*cr.getSpacing();
      //      _capacity += sites*crSite.height*crSite.width;
      minSiteHeight = min(minSiteHeight, crSite.height);
      minSiteSpacing = min(minSiteSpacing, cr.getSpacing());
      avgSiteHeight += sites * crSite.height;
      avgSiteSpacing += sites * cr.getSpacing();

      double leftEdge = cr[_startOffsets[r].first].getXStart();
      leftEdge += _startOffsets[r].second * cr.getSpacing();
      double rightEdge = cr[_endOffsets[r].first].getXStart();
      rightEdge += _endOffsets[r].second * cr.getSpacing() + cr.getSiteWidth();

      _bbox += Point(leftEdge, cr.getYCoord());
      _bbox += Point(rightEdge, cr.getYCoord() + cr.getHeight());
    }
    _capacity *= minSiteHeight;
    // avgSiteHeight /=_numSites;
    // avgSiteSpacing/=_numSites;
  }
  if (_capo.getParams().cogLocation) computeCenterOfGravity();

  _maxRecCellArea = min(_maxRecCellArea, _capacity);
}

void CapoBin::computeCenterOfGravity(void) {
  // cout<<endl<<"DPDEBUG computing center of gravity"<<endl;
  // cout<<"Before:"<<endl;
  // cout<<"\tBBox:"<<_bbox<<endl;
  // cout<<"\tCOG:"<<_centerOfGravity<<endl;

  // Does this algorithm properly account for subrows?
  _bbox = BBox();  // reset BBox to "empty"
  double xsum = 0.0, ysum = 0.0;
  unsigned Nsites = 0;
  for (unsigned r = 0; r != _coreRows.size(); r++) {
    const RBPCoreRow& cr = *_coreRows[r];
    // const RBPSite& crSite = cr.getSite();

    unsigned sites = getContainedSitesInCoreRow(r);
    Nsites += sites;

    double leftEdge = cr[_startOffsets[r].first].getXStart();
    leftEdge += _startOffsets[r].second * cr.getSpacing();
    double rightEdge = cr[_endOffsets[r].first].getXStart();
    rightEdge += _endOffsets[r].second * cr.getSpacing() + cr.getSiteWidth();

    _bbox += Point(leftEdge, cr.getYCoord());
    _bbox += Point(rightEdge, cr.getYCoord() + cr.getHeight());
    xsum += sites * (0.5 * fabs(rightEdge - leftEdge) + leftEdge);
    ysum += sites * cr.getYCoord();
  }

  // assigning the center of gravity to the geom center to test
  // just bbox triming
  //_centerOfGravity = Point( xsum / _numSites, ysum / _numSites);
  _centerOfGravity = _bbox.getGeomCenter();
  // cout<<"After:"<<endl;
  // cout<<"\tBBox:"<<_bbox<<endl;
  // cout<<"\tCOG:"<<_centerOfGravity<<endl<<endl;
}

Point CapoBin::getCenter() const {
  if (_capo.getParams().cogLocation) {
    abkfatal(_centerOfGravity.x != DBL_MAX && _centerOfGravity.y != DBL_MAX,
             "getCenter() called but _centerOfGravity not initialized.");
    return _centerOfGravity;
  } else
    return _bbox.getGeomCenter();
}

bool CapoBin::shouldAlignCL2Obs(double& splitPt, bool& splitHoriz) {
  double binArea = _bbox.getArea();
  double layoutArea = _capo.getRBP().getBBox(false).getArea();
  BBox nodeBBox;
  BBox largestObs;
  BBox p0BBox = _bbox;
  bool foundLargest = false;
  double largestArea = -DBL_MAX;
  double p0Area, p1Area, minRelDiff, maxRelDiff;
  const PlacementWOrient& pl = _capo.getPlacement();

  if (binArea / layoutArea < 0.0625) return false;

  const HGraphWDimensions& hgraph = _capo.getNetlistHGraph();

  for (unsigned i = 0; i < _obstacleCellIds.size(); ++i) {
    unsigned index = _obstacleCellIds[i];
    unsigned angle = pl.getOrient(index).getAngle();
    bool notRotated = angle == 0 || angle == 180;
    double nodeWidth = notRotated ? (hgraph.getNodeWidth(index))
                                  : (hgraph.getNodeHeight(index));
    double nodeHeight = notRotated ? (hgraph.getNodeHeight(index))
                                   : (hgraph.getNodeWidth(index));
    nodeBBox.clear();
    nodeBBox.add(pl[index].x, pl[index].y);
    nodeBBox.add(pl[index].x + nodeWidth, pl[index].y + nodeHeight);

    BBox intersection = nodeBBox.intersect(_bbox);
    double intersectionArea = intersection.getArea();
    if ((intersectionArea / binArea) > 0.0625) {
      if (intersectionArea > largestArea) {
        foundLargest = true;
        largestArea = intersectionArea;
        largestObs = intersection;
      }
    }
  }

  if (foundLargest) {
    bool validObs = false;
    if ((largestObs.xMin > _bbox.xMin || largestObs.xMax < _bbox.xMax)) {
      if ((largestObs.getHeight() <= largestObs.getWidth())) {
        if ((largestObs.yMin > _bbox.yMin || largestObs.yMax < _bbox.yMax) &&
            (equalDouble(largestObs.xMin, _bbox.xMin) ||
             equalDouble(largestObs.xMax, _bbox.xMax)) &&
            (_bbox.getWidth() <= 4 * _bbox.getHeight())) {
          splitHoriz = true;
          validObs = true;
        } else if (_bbox.getHeight() <= 4 * _bbox.getWidth()) {
          splitHoriz = false;
          validObs = true;
        } else if ((largestObs.yMin > _bbox.yMin ||
                    largestObs.yMax < _bbox.yMax)) {
          splitHoriz = true;
          validObs = true;
        }
      }
    }

    if ((largestObs.yMin > _bbox.yMin || largestObs.yMax < _bbox.yMax)) {
      if ((largestObs.getHeight() >= largestObs.getWidth())) {
        if ((largestObs.xMin > _bbox.xMin || largestObs.xMax < _bbox.xMax) &&
            (equalDouble(largestObs.yMin, _bbox.yMin) ||
             equalDouble(largestObs.yMax, _bbox.yMax)) &&
            (_bbox.getHeight() <= 4 * _bbox.getWidth())) {
          splitHoriz = false;
          validObs = true;
        } else if (_bbox.getWidth() <= 4 * _bbox.getHeight()) {
          splitHoriz = true;
          validObs = true;
        } else if ((largestObs.xMin > _bbox.xMin ||
                    largestObs.xMax < _bbox.xMax)) {
          splitHoriz = false;
          validObs = true;
        }
      }
    }

    if (validObs) {
      if (splitHoriz) {
        p0BBox.yMin = largestObs.yMin;
        p0Area = getContainedAreaInBBox(p0BBox);
        p1Area = _capacity - p0Area;
        minRelDiff = fabs(p0Area - p1Area) / _capacity;
        p0BBox.yMin = _bbox.yMin;

        p0BBox.yMin = largestObs.yMax;
        p0Area = getContainedAreaInBBox(p0BBox);
        p1Area = _capacity - p0Area;
        maxRelDiff = fabs(p0Area - p1Area) / _capacity;
        p0BBox.yMin = _bbox.yMin;

        /*Objective here is different: to have the maximum
          relative difference in whitespace to isolate slivers
          early on*/
        if (minRelDiff > maxRelDiff && largestObs.yMin > _bbox.yMin)
          splitPt = largestObs.yMin;
        else if (largestObs.yMax < _bbox.yMax)
          splitPt = largestObs.yMax;
        else if (largestObs.yMin > _bbox.yMin)
          splitPt = largestObs.yMin;
        else
          foundLargest = false;

        if (foundLargest) splitPt = 0.5 * (largestObs.yMax + largestObs.yMin);
      } else {
        p0BBox.xMin = largestObs.xMin;
        p0Area = getContainedAreaInBBox(p0BBox);
        p1Area = _capacity - p0Area;
        minRelDiff = fabs(p0Area - p1Area) / _capacity;
        p0BBox.xMin = _bbox.xMin;

        p0BBox.xMin = largestObs.xMax;
        p0Area = getContainedAreaInBBox(p0BBox);
        p1Area = _capacity - p0Area;
        maxRelDiff = fabs(p0Area - p1Area) / _capacity;
        p0BBox.xMin = _bbox.xMin;

        /*Objective here is different: to have the maximum
          relative difference in whitespace to isolate slivers
          early on*/
        if (minRelDiff > maxRelDiff && largestObs.xMin > _bbox.xMin)
          splitPt = largestObs.xMin;
        else if (largestObs.xMax < _bbox.xMax)
          splitPt = largestObs.xMax;
        else if (largestObs.xMin > _bbox.xMin)
          splitPt = largestObs.xMin;
        else
          foundLargest = false;

        if (foundLargest) splitPt = 0.5 * (largestObs.xMax + largestObs.xMin);
      }
    } else {
      foundLargest = false;
    }
  }
  if (foundLargest) {
    cout << "Bin : " << _bbox << " Obstruction: " << largestObs << endl;
    cout << "Trying to split the bin ";
    if (splitHoriz)
      cout << "horizontally";
    else
      cout << "vertically";
    cout << " at location " << splitPt << endl;
  } else {
    // cout<<"Found no big obstacle \n";
  }
  return (foundLargest);
}

void CapoBin::setIndexToNextIndex() { _index = _capo.nextBinNum++; }

void CapoBin::decrementIdx() { --_capo.nextBinNum; }

CapoBin::~CapoBin() {
  unLinkNeighbors();
  //  if (_savedSoln != NULL) delete _savedSoln; _savedSoln = NULL;
}

CapoBin::hCutLineInfo CapoBin::scanCutLinesH(const SubHGraph& probSubHGraph,
                                             double p0Min, double p0Max,
                                             double p1Min, double p1Max) const {
  unsigned bestSplitRow = UINT_MAX;
  double bestCut = DBL_MAX;
  double bestIllegality = DBL_MAX;
  double bestImbalance = DBL_MAX;
  double bestP0Area = DBL_MAX;
  double bestP1Area = DBL_MAX;

  const PlacementWOrient& oraclePl = _capo.getOraclePlacement();
  const HGraphWDimensions& nlhgraph = _capo.getNetlistHGraph();

  vector<pair<double, double> > cellYCenter(_cellIds.size());

  for (unsigned i = 0; i < _cellIds.size(); ++i) {
    unsigned angle = oraclePl.getOrient(_cellIds[i]).getAngle();
    bool notRotated = angle == 0 || angle == 180;
    double cellHeight = notRotated ? nlhgraph.getNodeHeight(_cellIds[i])
                                   : nlhgraph.getNodeWidth(_cellIds[i]);

    cellYCenter[i].first = oraclePl[_cellIds[i]].y + 0.5 * cellHeight;
    cellYCenter[i].second = nlhgraph.getWeight(_cellIds[i]);
  }

  sort(cellYCenter.begin(), cellYCenter.end());

  // for nets, need the min and max y coords of the bbox of the net
  vector<pair<double, double> > netMinY(probSubHGraph.getNumEdges());
  vector<pair<double, double> > netMaxY(probSubHGraph.getNumEdges());

  for (unsigned i = 0; i < probSubHGraph.getNumEdges(); ++i) {
    const HGFEdge& edge = probSubHGraph.getEdgeByIdx(i);

    netMinY[i].first = DBL_MAX;
    netMaxY[i].first = -DBL_MAX;

    for (itHGFNodeLocal n = edge.nodesBegin(); n != edge.nodesEnd(); ++n) {
      unsigned newNodeIdx = (**n).getIndex();
      if (newNodeIdx == 0) {
        netMaxY[i].first = DBL_MAX;
      } else if (newNodeIdx == 1) {
        netMinY[i].first = -DBL_MAX;
      } else {
        unsigned origNodeIdx = probSubHGraph.newNode2OrigIdx(newNodeIdx);

        unsigned angle = oraclePl.getOrient(origNodeIdx).getAngle();
        bool notRotated = angle == 0 || angle == 180;
        double cellHeight = notRotated ? nlhgraph.getNodeHeight(origNodeIdx)
                                       : nlhgraph.getNodeWidth(origNodeIdx);

        netMinY[i].first =
            min(netMinY[i].first, oraclePl[origNodeIdx].y + 0.5 * cellHeight);
        netMaxY[i].first =
            max(netMaxY[i].first, oraclePl[origNodeIdx].y + 0.5 * cellHeight);
      }
    }

    // the net weight is rounded because FMPart rounds all of its edge weights
    netMinY[i].second = netMaxY[i].second = rint(edge.getWeight());
  }

  sort(netMinY.begin(), netMinY.end());
  sort(netMaxY.begin(), netMaxY.end());

  unsigned cellCenterIdx = 0;
  unsigned netMinIdx = 0, netMaxIdx = 0;

  double curNetCut = 0.;
  double p0Area = getTotalCellArea();
  double p1Area = 0.;
  for (unsigned splitRow = 0; splitRow < _coreRows.size(); ++splitRow) {
    double splitRowCoord = _coreRows[splitRow]->getYCoord();

    while (cellCenterIdx < cellYCenter.size() &&
           lessOrEqualDouble(cellYCenter[cellCenterIdx].first, splitRowCoord)) {
      p0Area -= cellYCenter[cellCenterIdx].second;
      p1Area += cellYCenter[cellCenterIdx].second;
      ++cellCenterIdx;
    }

    while (netMinIdx < netMinY.size() &&
           lessThanDouble(netMinY[netMinIdx].first, splitRowCoord)) {
      curNetCut += netMinY[netMinIdx].second;
      ++netMinIdx;
    }

    while (netMaxIdx < netMaxY.size() &&
           lessOrEqualDouble(netMaxY[netMaxIdx].first, splitRowCoord)) {
      curNetCut -= netMaxY[netMaxIdx].second;
      ++netMaxIdx;
    }

    if ((splitRow >= (getNumRows() + 1) / 3) &&
        (splitRow <= (2 * getNumRows() - 1) / 3)) {
      double imbalance = fabs(p0Area - p1Area);
      if (p0Min <= p0Area && p0Area <= p0Max && p1Min <= p1Area &&
          p1Area <= p1Max) {
        if (curNetCut < bestCut ||
            (curNetCut == bestCut && imbalance < bestImbalance)) {
          bestCut = curNetCut;
          bestSplitRow = splitRow;
          bestIllegality = 0.;
          bestImbalance = imbalance;
          bestP0Area = p0Area;
          bestP1Area = p1Area;
        }
      } else if (bestIllegality > 0) {
        double p0Illegality = max(p0Min - p0Area, p0Area - p0Max);
        double p1Illegality = max(p1Min - p1Area, p1Area - p1Max);
        double illegality = max(p0Illegality, p1Illegality);
        if (illegality < bestIllegality) {
          bestCut = curNetCut;
          bestSplitRow = splitRow;
          bestIllegality = illegality;
          bestImbalance = imbalance;
          bestP0Area = p0Area;
          bestP1Area = p1Area;
        }
      }
    }
  }

  CapoBin::hCutLineInfo bestCutLine;
  bestCutLine.splitRow = bestSplitRow;
  bestCutLine.cutCost = bestCut;
  bestCutLine.p0Area = bestP0Area;
  bestCutLine.p1Area = bestP1Area;

  return bestCutLine;
}

CapoBin::vCutLineInfo CapoBin::scanCutLinesV(const SubHGraph& probSubHGraph,
                                             double p0Min, double p0Max,
                                             double p1Min, double p1Max) const {
  double bestXSplit = DBL_MAX;
  double bestCut = DBL_MAX;
  double bestIllegality = DBL_MAX;
  double bestImbalance = DBL_MAX;
  double bestP0Area = DBL_MAX;
  double bestP1Area = DBL_MAX;

  const PlacementWOrient& oraclePl = _capo.getOraclePlacement();
  const HGraphWDimensions& nlhgraph = _capo.getNetlistHGraph();

  vector<pair<double, double> > cellXCenter(_cellIds.size());

  for (unsigned i = 0; i < _cellIds.size(); ++i) {
    unsigned angle = oraclePl.getOrient(_cellIds[i]).getAngle();
    bool notRotated = angle == 0 || angle == 180;
    double cellWidth = notRotated ? nlhgraph.getNodeWidth(_cellIds[i])
                                  : nlhgraph.getNodeHeight(_cellIds[i]);

    cellXCenter[i].first = oraclePl[_cellIds[i]].x + 0.5 * cellWidth;
    cellXCenter[i].second = nlhgraph.getWeight(_cellIds[i]);
  }

  sort(cellXCenter.begin(), cellXCenter.end());

  // for nets, need the min and max y coords of the bbox of the net
  vector<pair<double, double> > netMinX(probSubHGraph.getNumEdges());
  vector<pair<double, double> > netMaxX(probSubHGraph.getNumEdges());

  for (unsigned i = 0; i < probSubHGraph.getNumEdges(); ++i) {
    const HGFEdge& edge = probSubHGraph.getEdgeByIdx(i);

    netMinX[i].first = DBL_MAX;
    netMaxX[i].first = -DBL_MAX;

    for (itHGFNodeLocal n = edge.nodesBegin(); n != edge.nodesEnd(); ++n) {
      unsigned newNodeIdx = (**n).getIndex();
      if (newNodeIdx == 0) {
        netMinX[i].first = -DBL_MAX;
      } else if (newNodeIdx == 1) {
        netMaxX[i].first = DBL_MAX;
      } else {
        unsigned origNodeIdx = probSubHGraph.newNode2OrigIdx(newNodeIdx);

        unsigned angle = oraclePl.getOrient(origNodeIdx).getAngle();
        bool notRotated = angle == 0 || angle == 180;
        double cellWidth = notRotated ? nlhgraph.getNodeWidth(origNodeIdx)
                                      : nlhgraph.getNodeHeight(origNodeIdx);

        netMinX[i].first =
            min(netMinX[i].first, oraclePl[origNodeIdx].x + 0.5 * cellWidth);
        netMaxX[i].first =
            max(netMaxX[i].first, oraclePl[origNodeIdx].x + 0.5 * cellWidth);
      }
    }

    // the net weight is rounded because FMPart rounds all of its edge weights
    netMinX[i].second = netMaxX[i].second = rint(edge.getWeight());
  }

  sort(netMinX.begin(), netMinX.end());
  sort(netMaxX.begin(), netMaxX.end());

  unsigned cellCenterIdx = 0;
  unsigned netMinIdx = 0, netMaxIdx = 0;

  double curNetCut = 0.;
  double p0Area = 0.;
  double p1Area = getTotalCellArea();
  double spacing = _coreRows[0]->getSpacing();
  vector<pair<double, double> > areaBalances;

  for (double xSplit = _bbox.xMin; xSplit <= _bbox.xMax; xSplit += spacing) {
    while (cellCenterIdx < cellXCenter.size() &&
           lessOrEqualDouble(cellXCenter[cellCenterIdx].first, xSplit)) {
      p0Area += cellXCenter[cellCenterIdx].second;
      p1Area -= cellXCenter[cellCenterIdx].second;
      ++cellCenterIdx;
    }

    while (netMinIdx < netMinX.size() &&
           lessOrEqualDouble(netMinX[netMinIdx].first, xSplit)) {
      curNetCut += netMinX[netMinIdx].second;
      ++netMinIdx;
    }

    while (netMaxIdx < netMaxX.size() &&
           lessOrEqualDouble(netMaxX[netMaxIdx].first, xSplit)) {
      curNetCut -= netMaxX[netMaxIdx].second;
      ++netMaxIdx;
    }

    if ((xSplit >= _bbox.xMin + 0.3333 * (_bbox.xMax - _bbox.xMin)) &&
        (xSplit <= _bbox.xMin + 0.6667 * (_bbox.xMax - _bbox.xMin))) {
      double imbalance = fabs(p0Area - p1Area);
      if (p0Min <= p0Area && p0Area <= p0Max && p1Min <= p1Area &&
          p1Area <= p1Max) {
        if (curNetCut < bestCut ||
            (curNetCut == bestCut && imbalance < bestImbalance)) {
          bestCut = curNetCut;
          bestXSplit = xSplit;
          bestIllegality = 0.;
          bestImbalance = imbalance;
          bestP0Area = p0Area;
          bestP1Area = p1Area;
        }
      } else if (bestIllegality > 0) {
        double p0Illegality = max(p0Min - p0Area, p0Area - p0Max);
        double p1Illegality = max(p1Min - p1Area, p1Area - p1Max);
        double illegality = max(p0Illegality, p1Illegality);
        if (illegality < bestIllegality) {
          bestCut = curNetCut;
          bestXSplit = xSplit;
          bestIllegality = illegality;
          bestImbalance = imbalance;
          bestP0Area = p0Area;
          bestP1Area = p1Area;
        }
      }
    }
  }

  CapoBin::vCutLineInfo bestCutLine;
  bestCutLine.xSplit = bestXSplit;
  bestCutLine.cutCost = bestCut;
  bestCutLine.p0Area = bestP0Area;
  bestCutLine.p1Area = bestP1Area;

  return bestCutLine;
}
