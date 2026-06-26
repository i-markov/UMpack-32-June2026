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

// created by Andrew Caldwell on 02/02/99
#include "hgWDims.h"

#include <cmath>
#include <string>

#include "ABKCommon/abkassert.h"
#include "ABKCommon/pathDelims.h"

using std::cerr;
using std::cout;
using std::endl;
using std::vector;
using std::string;

const float HGraphWDimensions::DEFAULT_MAX_AR = 2.0;
const float HGraphWDimensions::DEFAULT_MIN_AR = 0.5;

HGraphWDimensions::HGraphWDimensions(unsigned numNodes, unsigned numWeights,
                                     HGraphParameters param)
    : SubHGraph(numNodes, numWeights, param) {
  setupUnitDimensions();
  setupCenterPinOffsets();
  setupEmptySymmetries();
  setupDefaultSoftNodes();
  _isMacro = vector<bool>(numNodes, false);
}

HGraphWDimensions::HGraphWDimensions(unsigned numNodes, unsigned numWeights,
                                     HGraphParameters param,
                                     const vector<float>& heights,
                                     const vector<float>& widths)
    : SubHGraph(numNodes, numWeights, param),
      _nodeHeights(heights),
      _nodeWidths(widths) {
  setupUnitDimensions();
  setupCenterPinOffsets();
  setupEmptySymmetries();
  setupDefaultSoftNodes();
  _isMacro = vector<bool>(numNodes, false);
}

HGraphWDimensions::HGraphWDimensions(const HGraphWDimensions& orig)
    : SubHGraph(orig),
      _nodeHeights(orig._nodeHeights),
      _nodeWidths(orig._nodeWidths),
      _nodeSymmetries(orig._nodeSymmetries),
      _isMacro(orig._isMacro),
      _isSoft(orig._isSoft),
      _nodeMaxAR(orig._nodeMaxAR),
      _nodeMinAR(orig._nodeMinAR),
      _pinOffsets(orig._pinOffsets),
      _delayInfo(orig._delayInfo),
      _nodesOnEdgePinIdx(orig._nodesOnEdgePinIdx),
      _edgesOnNodePinIdx(orig._edgesOnNodePinIdx),
      _noeBeginEnd(orig._noeBeginEnd),
      _eonBeginEnd(orig._eonBeginEnd),
#ifdef SIGNAL_DIRECTIONS
      _srcsWPins(orig._srcsWPins),
      _snksWPins(orig._snksWPins),
#endif
      _srcSnksWPins(orig._srcSnksWPins) {
}

HGraphWDimensions::HGraphWDimensions(const char* fileName1,
                                     const char* fileName2,
                                     const char* fileName3,
                                     HGraphParameters param)
    //   :HGraphFixed(param),
    : SubHGraph(param) {
  abkfatal(fileName1, "null char* for first filename in HGraphWDims ctor");

  const char* auxFileName = 0;
  const char* netDFileName = 0;
  const char* areMFileName = 0;
  const char* wtsFileName = 0;
  const char* netsFileName = 0;
  const char* nodesFileName = 0;

  const char* files[3] = {fileName1, fileName2, fileName3};

  unsigned f;
  for (f = 0; f < 3; f++) {
    if (!files[f]) continue;

    if (strstr(files[f], ".aux"))
      auxFileName = files[f];
    else if (strstr(files[f], ".wts"))
      wtsFileName = files[f];
    else if (strstr(files[f], ".are"))
      areMFileName = files[f];
    else if (strstr(files[f], ".nodes"))
      nodesFileName = files[f];
    else if (strstr(files[f], ".netD"))
      netDFileName = files[f];
    else if (strstr(files[f], ".nets"))
      netsFileName = files[f];
    else if (strstr(files[f], ".net"))
      netDFileName = files[f];
  }

  abkfatal(!((netDFileName || areMFileName) &&
             (netsFileName || nodesFileName || wtsFileName)),
           "can't mix netD/areM and nets/nodes/wts");

  if (auxFileName)
    parseAux(auxFileName);
  else if (netDFileName) {
    readNetD(netDFileName);
    if (areMFileName)
      readAreM(areMFileName);
    else {  // set unit weights
      abkfatal(_numMultiWeights == 0, "there should be no weights here");
      _numMultiWeights = 1;
      _multiWeights = std::vector<HGWeight>(_nodes.size(), 1.0);
    }
  } else if (nodesFileName) {
    abkfatal(netsFileName, "must have both a nodes and nets file");

    readNodes(nodesFileName);
    readNets(netsFileName);
    if (wtsFileName)
      readWts(wtsFileName);
    else {  // set unit weights
      abkfatal(_numMultiWeights == 0, "there should be no weights here");
      _numMultiWeights = 1;
      _multiWeights = std::vector<HGWeight>(_nodes.size(), 1.0);
    }
  } else {
    cerr << "aux, netD or nets file not found.  ";
    cerr << "Filename arguments were:" << endl;
    if (fileName1) cout << fileName1 << endl;
    if (fileName2) cout << fileName2 << endl;
    if (fileName3) cout << fileName3 << endl;
  }

  finalize();

  setupUnitDimensions();
  setupCenterPinOffsets();
  setupEmptySymmetries();
  setupDefaultSoftNodes();
  _isMacro = vector<bool>(_nodes.size(), false);
}

HGraphWDimensions::HGraphWDimensions(const HGraphFixed& hg,
                                     const vector<float>& heights,
                                     const vector<float>& widths)

    //   :HGraphFixed(hg),
    : SubHGraph(hg), _nodeHeights(heights), _nodeWidths(widths) {
  setupUnitDimensions();
  setupCenterPinOffsets();
  setupEmptySymmetries();
  setupDefaultSoftNodes();
  _isMacro = vector<bool>(_nodes.size(), false);
}

void HGraphWDimensions::setupUnitDimensions() {
  if (_nodeWidths.size() == 0) {
    _nodeWidths = vector<float>(_nodes.size(), 1.0);
    _nodeHeights = vector<float>(_nodes.size(), 1.0);
  }
}

void HGraphWDimensions::setupEmptySymmetries() {
  if (_nodeSymmetries.size() != _nodes.size())
    _nodeSymmetries = vector<Symmetry>(_nodes.size());
}

void HGraphWDimensions::setupCenterPinOffsets() {
  if (_pinOffsets.size() == 0) {
    _nodesOnEdgePinIdx = vector<unsigned>(_numPins);
    _edgesOnNodePinIdx = vector<unsigned>(_numPins);
    _noeBeginEnd = vector<unsigned>(_edges.size() + 1);
    _eonBeginEnd = vector<unsigned>(_nodes.size() + 1);
    _pinOffsets = vector<PinPoint>(_nodes.size());

    unsigned offset = 0;
    unsigned n;
    for (n = 0; n < _nodes.size(); n++) {
      _eonBeginEnd[n] = offset;
      _pinOffsets[n].offset.x = _nodeWidths[n] / 2.0;
      _pinOffsets[n].offset.y = _nodeHeights[n] / 2.0;

      for (unsigned e = 0; e < _nodes[n]->getDegree(); e++, offset++)
        _edgesOnNodePinIdx[offset] = n;  // all pins are in the middle
    }
    _eonBeginEnd[n] = offset;

    itHGFNodeLocal nIt;
    offset = 0;
    unsigned e;
    for (e = 0; e < _edges.size(); e++) {
      _noeBeginEnd[e] = offset;
      for (nIt = _edges[e]->nodesBegin(); nIt != _edges[e]->nodesEnd();
           nIt++, offset++)
        _nodesOnEdgePinIdx[offset] = (*nIt)->getIndex();
    }
    _noeBeginEnd[e] = offset;
  }
}

void HGraphWDimensions::makeCenterPinOffsets() {
  unsigned n;
  for (n = 0; n < _nodes.size(); n++) {
    for (unsigned e = 0; e < _nodes[n]->getDegree(); e++) {
      unsigned pOffsetIdx = _edgesOnNodePinIdx[_eonBeginEnd[n] + e];
      _pinOffsets[pOffsetIdx].offset.x = _nodeWidths[n] / 2.0;
      _pinOffsets[pOffsetIdx].offset.y = _nodeHeights[n] / 2.0;
    }
  }
}

void HGraphWDimensions::updateNOEPinOffsets(unsigned nodeOffset,
                                            unsigned edgeId,
                                            Point newPinOffset) {
  unsigned pOffsetIdx = _nodesOnEdgePinIdx[_noeBeginEnd[edgeId] + nodeOffset];
  _pinOffsets[pOffsetIdx].offset = newPinOffset;
}

void HGraphWDimensions::updateEONPinOffsets(unsigned edgeOffset,
                                            unsigned nodeId,
                                            Point newPinOffset) {
  unsigned pOffsetIdx = _edgesOnNodePinIdx[_eonBeginEnd[nodeId] + edgeOffset];
  _pinOffsets[pOffsetIdx].offset = newPinOffset;
}

void HGraphWDimensions::calcRealOffset(unsigned cellId, Point& pinOffset,
                                       const Orient& ori) const {
  float h = getNodeHeight(cellId);
  float w = getNodeWidth(cellId);

  if (!ori.isFaceUp()) pinOffset.x = w - pinOffset.x;

  float tmp;
  switch (ori.getAngle()) {
    case 0:  //( x,  y)
      break;
    case 90:  //( y, -x)
      tmp = pinOffset.x;
      pinOffset.x = pinOffset.y;
      pinOffset.y = w - tmp;
      break;
    case 180:  //(-x, -y)
      pinOffset.x = w - pinOffset.x;
      pinOffset.y = h - pinOffset.y;
      break;
    case 270:  //(-y,  x)
      tmp = pinOffset.y;
      pinOffset.y = pinOffset.x;
      pinOffset.x = h - tmp;
      break;
    default:
      abkfatal(0, " Orientation corrupted ");
  }
}

bool HGraphWDimensions::isConsistent() const {
  abkfatal(_finalized, "HGraphWDims was not finalized");
  return HGraphFixed::isConsistent();
}

void HGraphWDimensions::finalize() {
  if (_param.verb.getForActions() > 5) cout << "Finalizing HGWDims" << endl;
  abkfatal(_finalized == false, "cannot run finalize twice");

#ifdef SIGNAL_DIRECTIONS
  if ((!_srcs.empty() || !_snks.empty() || !_srcSnks.empty()) &&
      (_srcsWPins.size() || _snksWPins.size() || _srcSnksWPins.size()))
#else
  if (!_srcSnks.empty() && _srcSnksWPins.size())
#endif
    abkfatal(0, "Error: both HGraph and HGWDims node-net pairs are populated");

#ifdef SIGNAL_DIRECTIONS
  if (!_srcs.empty() || !_snks.empty() || !_srcSnks.empty())
#else
  if (!_srcSnks.empty())
#endif
  {
    HGraphFixed::finalize();
    createIdentityMappings();
    return;
  }

  if (_addEdgeStyle == FAST_ADDEDGE_USED) {
#ifdef SIGNAL_DIRECTIONS
    for (unsigned v = 0; v != getNumNodes(); ++v) {
      HGFNode& node = getNodeByIdx(v);
      node._srcSnksBegin = node._edges.begin();
      node._snksBegin = node._edges.end();
    }
    for (unsigned e = 0; e != getNumEdges(); ++e) {
      HGFEdge& edge = getEdgeByIdx(e);
      edge._srcSnksBegin = edge._nodes.begin();
      edge._snksBegin = edge._nodes.end();
    }
#endif
    _finalized = true;
    createIdentityMappings();
    return;
  }

  /*    abkwarn(!_param.removeBigNets,
        "HGraphWDims does not support removal of large edges");*/

  vector<unsigned> nodeDegrees(getNumNodes());
  vector<unsigned> edgeDegrees(getNumEdges());

  vector<HGWDimsPin>::iterator it;
  _numPins = 0;

#ifdef SIGNAL_DIRECTIONS
  for (it = _srcsWPins.begin(); it != _srcsWPins.end(); ++it) {
    _numPins++;
    ++nodeDegrees[it->nodeId];
    ++edgeDegrees[it->edgeId];
  }
#endif

  for (it = _srcSnksWPins.begin(); it != _srcSnksWPins.end(); ++it) {
    _numPins++;
    ++nodeDegrees[it->nodeId];
    ++edgeDegrees[it->edgeId];
  }

#ifdef SIGNAL_DIRECTIONS
  for (it = _snksWPins.begin(); it != _snksWPins.end(); ++it) {
    _numPins++;
    ++nodeDegrees[it->nodeId];
    ++edgeDegrees[it->edgeId];
  }
#endif

  if (_param.removeBigNets) {
    cout << "Removing nets..." << endl;
    // std::list<pair<unsigned, unsigned> >::iterator it;
    vector<unsigned> old2new(getNumEdges());
    unsigned last = 0;
    unsigned old;
    unsigned rmCt = 0;
    for (old = 0; old != getNumEdges(); ++old) {
      if (edgeDegrees[old] <= _param.netThreshold)
        old2new[old] = last++;
      else {
        old2new[old] = UINT_MAX;
        rmCt++;
      }
    }

#ifdef SIGNAL_DIRECTIONS
    for (it = _srcsWPins.begin(); it != _srcsWPins.end(); ++it)
      (*it).edgeId = old2new[(*it).edgeId];
#endif

    for (it = _srcSnksWPins.begin(); it != _srcSnksWPins.end(); ++it)
      (*it).edgeId = old2new[(*it).edgeId];

#ifdef SIGNAL_DIRECTIONS
    for (it = _snksWPins.begin(); it != _snksWPins.end(); ++it)
      (*it).edgeId = old2new[(*it).edgeId];
#endif

    last = 0;
    for (unsigned i = 0; i != old2new.size(); ++i) {
      if (old2new[i] == UINT_MAX) {
        delete _edges[last];
        _edges.erase(_edges.begin() + last);
        //                _netNames.erase(_netNames.begin()+last);
        getNetNames().erase(getNetNames().begin() + last);
      } else {
        edgeDegrees[old2new[i]] = edgeDegrees[i];
        _edges[old2new[i]]->_index = old2new[i];
        last++;
      }
    }
    cout << "Thresholding " << rmCt << " Net(s) away." << endl;
  }

  // setup all pin offsets info, except the actual offsets.
  _pinOffsets.reserve(_numPins);
  _noeBeginEnd = vector<unsigned>(_edges.size() + 1);
  _eonBeginEnd = vector<unsigned>(_nodes.size() + 1);

  unsigned n, e;
  unsigned pinsBegin = 0;
  for (n = 0; n < _nodes.size(); n++) {
    _eonBeginEnd[n] = pinsBegin;
    pinsBegin += nodeDegrees[n];
    getNodeByIdx(n)._edges.reserve(nodeDegrees[n]);
  }
  _eonBeginEnd[n] = pinsBegin;

  pinsBegin = 0;
  for (e = 0; e < _edges.size(); e++) {
    _noeBeginEnd[e] = pinsBegin;
    pinsBegin += edgeDegrees[e];
    getEdgeByIdx(e)._nodes.reserve(edgeDegrees[e]);
  }
  _noeBeginEnd[e] = pinsBegin;

  _nodesOnEdgePinIdx = vector<unsigned>(_numPins);
  _edgesOnNodePinIdx = vector<unsigned>(_numPins);

  // these are used to put the pin offsets in the right loc
  vector<unsigned> nodesNextPin(_nodes.size(), 0);
  vector<unsigned> edgesNextPin(_edges.size(), 0);

#ifdef SIGNAL_DIRECTIONS
  for (it = _srcsWPins.begin(); it != _srcsWPins.end(); ++it) {
    if ((*it).edgeId == UINT_MAX) continue;
    if (edgeDegrees[(*it).edgeId] > _param.netThreshold) continue;
    unsigned nId = it->nodeId;
    unsigned eId = it->edgeId;

    getNodeByIdx(nId)._edges.push_back(&getEdgeByIdx(eId));
    getEdgeByIdx(eId)._nodes.push_back(&getNodeByIdx(nId));
    _nodesOnEdgePinIdx[_noeBeginEnd[eId] + edgesNextPin[eId]] =
        _edgesOnNodePinIdx[_eonBeginEnd[nId] + nodesNextPin[nId]] =
            _pinOffsets.size();
    edgesNextPin[eId]++;
    nodesNextPin[nId]++;
    _pinOffsets.push_back(it->pOffset);
  }
#endif

#ifdef SIGNAL_DIRECTIONS
  unsigned v;
  for (v = 0; v != getNumNodes(); ++v)
    getNodeByIdx(v)._srcSnksBegin = getNodeByIdx(v)._edges.end();
  for (e = 0; e != getNumEdges(); ++e)
    getEdgeByIdx(e)._srcSnksBegin = getEdgeByIdx(e)._nodes.end();
#endif

  for (it = _srcSnksWPins.begin(); it != _srcSnksWPins.end(); ++it) {
    if ((*it).edgeId == UINT_MAX) continue;
    if (edgeDegrees[(*it).edgeId] > _param.netThreshold) continue;
    unsigned nId = it->nodeId;
    unsigned eId = it->edgeId;
    getNodeByIdx(nId)._edges.push_back(&getEdgeByIdx(eId));
    getEdgeByIdx(eId)._nodes.push_back(&getNodeByIdx(nId));
    _nodesOnEdgePinIdx[_noeBeginEnd[eId] + edgesNextPin[eId]] =
        _edgesOnNodePinIdx[_eonBeginEnd[nId] + nodesNextPin[nId]] =
            _pinOffsets.size();
    edgesNextPin[eId]++;
    nodesNextPin[nId]++;
    _pinOffsets.push_back(it->pOffset);
  }

#ifdef SIGNAL_DIRECTIONS
  for (v = 0; v != getNumNodes(); ++v)
    getNodeByIdx(v)._snksBegin = getNodeByIdx(v)._edges.end();
  for (e = 0; e != getNumEdges(); ++e)
    getEdgeByIdx(e)._snksBegin = getEdgeByIdx(e)._nodes.end();
#endif

#ifdef SIGNAL_DIRECTIONS
  for (it = _snksWPins.begin(); it != _snksWPins.end(); ++it) {
    if ((*it).edgeId == UINT_MAX) continue;
    if (edgeDegrees[(*it).edgeId] > _param.netThreshold) continue;
    unsigned nId = it->nodeId;
    unsigned eId = it->edgeId;
    getNodeByIdx(nId)._edges.push_back(&getEdgeByIdx(eId));
    getEdgeByIdx(eId)._nodes.push_back(&getNodeByIdx(nId));
    _nodesOnEdgePinIdx[_noeBeginEnd[eId] + edgesNextPin[eId]] =
        _edgesOnNodePinIdx[_eonBeginEnd[nId] + nodesNextPin[nId]] =
            _pinOffsets.size();
    edgesNextPin[eId]++;
    nodesNextPin[nId]++;
    _pinOffsets.push_back(it->pOffset);
  }
#endif

#ifdef SIGNAL_DIRECTIONS
  _srcsWPins = vector<HGWDimsPin>(0);
  _snksWPins = vector<HGWDimsPin>(0);
#endif
  _srcSnksWPins = vector<HGWDimsPin>(0);

  sortNodes();
  sortEdges();

  _finalized = true;

  if (getNetNames().size() == 0) getNetNames().resize(_edges.size());

  createIdentityMappings();
}

void HGraphWDimensions::setNodeDims(float newHeight) {
  for (unsigned n = 0; n < _nodes.size(); n++) {
    if (getWeight(n) > 0) {
      _nodeHeights[n] = newHeight;
      _nodeWidths[n] = getWeight(n) / newHeight;
    }
  }
}

void HGraphWDimensions::setNodeDimsMacros(float newHeight) {
  float totArea = 0;
  unsigned n;
  for (n = 0; n < _nodes.size(); n++) totArea += getWeight(n);

  // float threshold = 0.001*totArea;
  for (n = 0; n < _nodes.size(); n++) {
    if (getWeight(n) > 0) {
      float potentialWidth = getWeight(n) / newHeight;

      // if(_nodes[n]->getWeight() < threshold)
      if (potentialWidth < 15 * newHeight) {
        _nodeHeights[n] = rint(newHeight);
        _nodeWidths[n] = rint(getWeight(n) / newHeight);
        _isMacro[n] = true;
      } else {
        float potentialHeight = rint(sqrt(getWeight(n)));
        float remainder = int(potentialHeight) % int(newHeight);
        if (remainder < newHeight / 2)
          potentialHeight -= remainder;
        else
          potentialHeight += (newHeight - remainder);

        _nodeHeights[n] = rint(potentialHeight);
        _nodeWidths[n] = rint(getWeight(n) / _nodeHeights[n]);
      }
    }
  }
}

void HGraphWDimensions::markNodesAsMacro(float maxWidth, float maxHeight) {
  for (unsigned i = 0; i < _nodes.size(); ++i) {
    if (!isTerminal(i)) {
      if (_nodeWidths[i] > maxWidth || _nodeHeights[i] > maxHeight)
        _isMacro[i] = true;
    }
  }
}

float HGraphWDimensions::getNodesArea() const {
  float nodesArea = 0;
  for (unsigned n = 0; n < _nodes.size(); n++) {
    if (!isTerminal(n)) nodesArea += _nodeWidths[n] * _nodeHeights[n];
  }
  return (nodesArea);
}

float HGraphWDimensions::getAvgNodeWidth(void) {
  float avgWidth = 0;
  int numNodes = 0;
  for (unsigned n = 0; n < _nodes.size(); n++) {
    if (!isTerminal(n)) {
      avgWidth += _nodeWidths[n];
      numNodes++;
    }
  }
  avgWidth /= numNodes;
  return (avgWidth);
}

float HGraphWDimensions::getAvgNodeHeight(void) {
  float avgHeight = 0;
  int numNodes = 0;
  for (unsigned n = 0; n < _nodes.size(); n++) {
    if (!isTerminal(n)) {
      avgHeight += _nodeHeights[n];
      ++numNodes;
    }
  }
  avgHeight /= numNodes;
  return (avgHeight);
}

void HGraphWDimensions::expandStdCellWidthToFitSiteWidth(double siteWidth) {
  float siteWidthF = static_cast<float>(siteWidth);

  for (unsigned i = getNumTerminals(); i < getNumNodes(); ++i) {
    if (!_isMacro[i]) {
      _nodeWidths[i] = ceilf(_nodeWidths[i] / siteWidthF) * siteWidthF;
    }
  }
}

void HGraphWDimensions::setNodeWeightsToAreas(double siteWidth,
                                              double siteHeight) {
  float siteWidthF = static_cast<float>(siteWidth);
  float siteHeightF = static_cast<float>(siteHeight);

  for (unsigned i = 0; i < getNumNodes(); ++i) {
    float width = ceilf(_nodeWidths[i] / siteWidthF) * siteWidthF;
    float height = ceilf(_nodeHeights[i] / siteHeightF) * siteHeightF;
    setWeight(i, width * height);
  }
}
