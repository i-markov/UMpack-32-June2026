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

// CHANGES
//  990511 aec	ver1.1	added getPinOffset functions that use an Orientation
//  990830 aec	ver2.0  updated to HGraph3.0 (parses nodes/nets,
//				supports multi-bin weights)
//  010721 ilm   ver3.1  added gate/wire delay info

#ifndef _HGRAPHWithDIMS_H_
#define _HGRAPHWithDIMS_H_

#include <cstddef>
#include <fstream>
#include <string>
#include <vector>

#include "ABKCommon/abkassert.h"
#include <string>
#include <vector>
#include "Geoms/bbox.h"
#include "Geoms/point.h"
#include "HGraph/hgFixed.h"
#include "HGraph/subHGraph.h"
#include "Placement/symOri.h"
#include "hgDelayInfo.h"

// HGraphWDimensions adds the following data to a HyperGraph
// 1) for each node, and height and width
// 2) for each node, the pin-offset for it's connection with adj nets
// 3) for each edge, the pin-offset and cellId's of adj nodes

class Pointf
/*
 * just like the point class, but with floats instead of doubles
 * XXX at some point, the Point class should be templatized
 */
{
 public:
  float x;
  float y;
  Pointf() : x(0.), y(0.) {}
  Pointf(float xx, float yy) : x(xx), y(yy) {}
  Pointf(const Pointf& orig) : x(orig.x), y(orig.y) {}
  Pointf(const Point& p)
      : x(static_cast<float>(p.x)), y(static_cast<float>(p.y)) {}
};

class PinPoint
// point (location) for a 'pin', which is derived from a set of related pins.
// a 'pin' is described by an (x,y) offset and the halfWidth and halfHeight
// of the bbox of the set of related pins.
{
 public:
  Pointf offset;
  Pointf halfBBoxDims;

  PinPoint() : offset(), halfBBoxDims() {}
  PinPoint(const Point& point) : offset(point), halfBBoxDims() {}
  PinPoint(const Pointf& offset_in, const Pointf& halfBBoxDims_in)
      : offset(offset_in), halfBBoxDims(halfBBoxDims_in) {}
  PinPoint(float xx, float yy, float hw, float hh)
      : offset(xx, yy), halfBBoxDims(hw, hh) {}
};

struct HGWDimsPin  // this stores nodeId, edgeId & pinOffset during
// construction
{
  HGWDimsPin() : nodeId(0), edgeId(0), pOffset() {}

  HGWDimsPin(unsigned n, unsigned e, const PinPoint& pt)
      : nodeId(n), edgeId(e), pOffset(pt) {}

  unsigned nodeId, edgeId;
  PinPoint pOffset;
};

// class HGraphWDimensions : public HGraphFixed
class HGraphWDimensions : public SubHGraph {
  friend class STA;

 protected:
  std::vector<float> _nodeHeights;        // height of each node
  std::vector<float> _nodeWidths;         // width  of each node
  std::vector<Symmetry> _nodeSymmetries;  // symmetry of each node
  std::vector<bool> _isMacro;             // is a node a macro

  // added by hhchan to keep track of soft block info
  std::vector<bool> _isSoft;      // is a block soft?
  std::vector<float> _nodeMaxAR;  // irrelevant unless _isSoft
  std::vector<float> _nodeMinAR;  // irrelevant unless _isSoft

  // collection of all pinOffsets, in no particular order
  std::vector<PinPoint> _pinOffsets;
  HGDelayInfo _delayInfo;  // includes a collection of Rs and Cs

  // the unsigneds are indices into _pinOffsets.
  // For each node, we store the beginning index for _edgesOnNodePinIdx,
  // in _eonBeginEnd,
  // and for each edge, the beginning index for _nodesOnEdgePinIdx
  // in _noeBeginEnd.
  // both have a sentinel value at the end (so that edgeId+1 gives
  // the id of the first pinoffset NOT in the edge)
  // the pin-offsets are given in terms of the offset of the nth
  // node on the edge.
  std::vector<unsigned> _nodesOnEdgePinIdx;
  std::vector<unsigned> _edgesOnNodePinIdx;

  std::vector<unsigned> _noeBeginEnd;
  std::vector<unsigned> _eonBeginEnd;

  // these are used in the expanded addS* and finalize functions
  // rather than the _src,_snk,_srcSnk lists the baseclass uses
#ifdef SIGNAL_DIRECTIONS
  std::vector<HGWDimsPin> _srcsWPins;
  std::vector<HGWDimsPin> _snksWPins;
#endif
  std::vector<HGWDimsPin> _srcSnksWPins;

  HGraphWDimensions(unsigned numNodes, unsigned numWeights,
                    HGraphParameters param = HGraphParameters());

  HGraphWDimensions(unsigned numNodes, unsigned numWeights,
                    HGraphParameters param, const std::vector<float>& heights,
                    const std::vector<float>& widths);

  // side-effect!! modifies pinOffset
  void calcRealOffset(unsigned nodeId, Point& pinOffset, const Orient&) const;

  virtual void readNets(const char* filename);
  virtual void readNodes(const char* filename);

  void setupUnitDimensions();
  void setupEmptySymmetries();
  void setupCenterPinOffsets();
  void setupDefaultSoftNodes();

#ifdef SIGNAL_DIRECTIONS
  void addSrcWPin(HGFNode& node, HGFEdge& edge, PinPoint pOffset);
  void addSnkWPin(HGFNode& node, HGFEdge& edge, PinPoint pOffset);
#endif
  void addSrcSnkWPin(HGFNode& node, HGFEdge& edge, PinPoint pOffset);

  virtual void finalize();

 public:
  HGraphWDimensions(const HGraphFixed& hg, const std::vector<float>& heights,
                    const std::vector<float>& widths);

  HGraphWDimensions(const HGraphWDimensions& orig);

  HGraphWDimensions(const char* fileName1, const char* fileName2 = NULL,
                    const char* fileName3 = NULL,
                    HGraphParameters param = HGraphParameters());

  virtual ~HGraphWDimensions() {}

  float getNodeHeight(unsigned id) const { return _nodeHeights[id]; }
  void setNodeWidth(float width, unsigned id) { _nodeWidths[id] = width; }
  void setNodeHeight(float height, unsigned id) { _nodeHeights[id] = height; }
  float getAvgNodeWidth(void);
  float getAvgNodeHeight(void);
  float getNodeWidth(unsigned id) const { return _nodeWidths[id]; }
  Symmetry getNodeSymmetry(unsigned id) const { return _nodeSymmetries[id]; }
  float getNodesArea() const;
  bool isNodeMacro(unsigned id) const { return _isMacro[id]; }
  void updateNodeMacroInfo(unsigned id, bool isMacro) {
    _isMacro[id] = isMacro;
  }
  // mark all nodes>maxArea as macros
  void markNodesAsMacro(float maxWidth, float maxHeight);

  // added by hhchan to keep track soft blocks
  static const float DEFAULT_MAX_AR;
  static const float DEFAULT_MIN_AR;

  float getNodeMaxAR(unsigned id) const { return _nodeMaxAR[id]; }
  float getNodeMinAR(unsigned id) const { return _nodeMinAR[id]; }
  bool isNodeSoft(unsigned id) const { return _isSoft[id]; }

  void setNodeMaxAR(float maxAR, unsigned id) { _nodeMaxAR[id] = maxAR; }
  void setNodeMinAR(float minAR, unsigned id) { _nodeMinAR[id] = minAR; }

  // side-effect: relocate the pins the center of their resp. nodes
  // warning: once a hard block -> soft, it's orig pin-offset is lost
  inline void updateNodeSoftInfo(unsigned id, float minAR = DEFAULT_MIN_AR,
                                 float maxAR = DEFAULT_MAX_AR);
  inline void shapeNodeProperly(unsigned id);
  inline void makeSoftNodesPinsCenter();
  inline void makeNodePinsCenter(unsigned id);
  inline void markMacrosAsSoft(float minAR = DEFAULT_MIN_AR,
                               float maxAR = DEFAULT_MAX_AR);
  void markSelectedAsSoftFromFile(const std::string& filename);

  const HGDelayInfo& getDelayInfo() const { return _delayInfo; }
  const std::vector<PinPoint>& getPinOffsets() const { return _pinOffsets; }
  bool haveDelayInfo() const { return _delayInfo.populated; }
  double getCIntC() const { return _delayInfo.cInt_c; }
  double getCIntK() const { return _delayInfo.cInt_k; }
  double getRIntC() const { return _delayInfo.rInt_c; }
  double getRIntK() const { return _delayInfo.rInt_k; }
  unsigned nodesOnEdgePinIdx(unsigned nodeOffset, unsigned edgeId) const;
  unsigned edgesOnNodePinIdx(unsigned edgeOffset, unsigned nodeId) const;
  Point getRealOffset(unsigned nodeId, const Point& pinOffset,
                      const Orient& ori) const;

  // sets all node's height to newHeight, and their widths to
  // their area(1st weight)/newHeight
  void setNodeDims(float newHeight);
  // sets all node's height to newHeight, and their widths to
  // their area(1st weight)/newHeight For large Macros use macros
  void setNodeDimsMacros(float newHeight);

  Point nodesOnEdgePinOffset(unsigned nOffset, unsigned eId,
                             Orient cellOri = Orient()) const;
  Point edgesOnNodePinOffset(unsigned eOffset, unsigned nId,
                             Orient cellOri = Orient()) const;
  BBox nodesOnEdgePinOffsetBBox(unsigned nOffset, unsigned eId,
                                Orient cellOri = Orient()) const;
  BBox edgesOnNodePinOffsetBBox(unsigned eOffset, unsigned nId,
                                Orient cellOri = Orient()) const;
  // double nodesOnEdgePinLoadCap(unsigned nOffset, unsigned eId) const;
  // double nodesOnEdgePinDrivRes(unsigned nOffset, unsigned eId) const;

  // double edgesOnNodePinLoadCap(unsigned eOffset, unsigned nId) const;
  // double edgesOnNodePinDrivRes(unsigned eOffset, unsigned nId) const;

  virtual bool saveAsNodesNetsWts(const char* baseFileName,
                                  bool fixMacros = false,
                                  bool saveNodeWts = true,
                                  bool saveNetWts = true) const;

  // by sadya to add dummy cells
  virtual void saveAsNodesNetsWtsWDummy(const char* baseFileName,
                                        float siteWidth, float siteHeight,
                                        float areaDummy = 0,
                                        bool fixMacros = false) const;
  // by sadya to add tethers to cells
  virtual void saveAsNodesNetsWtsWTether(
      const char* baseFileName, double fract,
      std::vector<unsigned>& tetheredCells, bool takeNodeConstrFrmFile = false,
      bool takeNetConstrFrmFile = false) const;
  void getCellsToTetherBFS(std::vector<unsigned>& cellToTether,
                           unsigned numTetherCells,
                           bool takeNodeConstrFrmFile = false,
                           bool takeNetConstrFrmFile = false) const;

  // added by sadya to shred the tall cells
  virtual void saveAsNodesNetsWtsShred(const char* baseFileName,
                                       float maxCoreRowHeight) const;
  virtual void saveAsNodesNetsWtsUnShred(const char* baseFileName,
                                         float maxCoreRowHeight) const;
  virtual void saveAsNodesNetsWtsShredHW(const char* baseFileName,
                                         float minWidthNode,
                                         float maxCoreRowHeight,
                                         unsigned singleNetWt = 0) const;

  // added by DP to functionalize the addition of nets
  // to allow for a flexible topology connecting shreds
  void AddNetsToShredsSingleNet(std::ofstream& nets, std::ofstream& wts,
                                float maxCoreRowHeight, float minWidthNode,
                                double maxWeight,
                                unsigned numNewSingleNets) const;
  void AddNetsToShredsGridTopology(std::ofstream& nets, std::ofstream& wts,
                                   float maxCoreRowHeight, float minWidthNode,
                                   double maxWeight) const;
  void AddNetsToShredsStarTopology(std::ofstream& nets, std::ofstream& wts,
                                   float maxCoreRowHeight, float minWidthNode,
                                   double maxWeight) const;

  virtual void saveAsNodesNetsFloorplan(const char* baseFileName) const;
  virtual void saveMacrosAsNodesNetsFloorplan(const char* baseFileName) const;

  // added by sadya to add nodes connected to congested regions
  virtual void saveAsNodesNetsWtsWCongInfo(
      const char* baseFileName, std::vector<double>& nodesXCongestion,
      std::vector<double>& nodesYCongestion, double whitespace,
      float rowHeight);

  virtual bool isConsistent() const;

  virtual void sortNodes();
  virtual void sortEdges();
  virtual void sortAsDB();

  // added by sadya to make existing poffsets at center
  void makeCenterPinOffsets();
  // added by sadya to update pinoffsets externally
  void updateNOEPinOffsets(unsigned nodeOffset, unsigned edgeId,
                           Point newPinOffset);
  void updateEONPinOffsets(unsigned edgeOffset, unsigned nodeId,
                           Point newPinOffset);

  // added by dp for easy manual construction
  // Instructions for manually constructing a HGraphWDimensions:
  // call createNewManualHGraphWDimensions to get an object to work with
  // for( all edges )
  //{
  //   Create an edge
  //   if(not all nodes on this edge are created)
  //     create missing nodes
  //   Add all nodes on the edge
  // }
  // done()
  static HGraphWDimensions* createNewManualHGraphWDimensions(unsigned);
  unsigned create_edge();
  unsigned create_edge(std::string name);
  unsigned create_node(std::string name, size_t width, size_t height,
                       bool term = false);
  unsigned create_node(size_t width, size_t height, bool term = false);
  void add_node_to_edge(unsigned node_idx, unsigned edge_idx, char direction);
  void done(void);

  void expandStdCellWidthToFitSiteWidth(double siteWidth);
  void setNodeWeightsToAreas(double siteWidth, double siteHeight);
};

inline unsigned HGraphWDimensions::nodesOnEdgePinIdx(unsigned nodeOffset,
                                                     unsigned edgeId) const {
  abkassert(edgeId < _edges.size(), "edgeId out of range in nodesOnEdgePinIdx");
  abkassert(nodeOffset < _edges[edgeId]->getDegree(),
            "nodeOffset of range in nodesOnEdgePinIdx");
  return _nodesOnEdgePinIdx[_noeBeginEnd[edgeId] + nodeOffset];
}

inline unsigned HGraphWDimensions::edgesOnNodePinIdx(unsigned edgeOffset,
                                                     unsigned nodeId) const {
  abkassert(nodeId < _nodes.size(),
            "nodeId out of range in edgesOnNodePinOffset");
  abkassert(edgeOffset < _nodes[nodeId]->getDegree(),
            "edgeOffset out of range in edgesOnNodePinOffset");

  return _edgesOnNodePinIdx[_eonBeginEnd[nodeId] + edgeOffset];
}

inline Point HGraphWDimensions::nodesOnEdgePinOffset(unsigned nodeOffset,
                                                     unsigned edgeId,
                                                     Orient cellOri) const {
  abkassert(edgeId < _edges.size(), "edgeId out of range in nodesOnEdgePinIdx");
  abkassert(nodeOffset < _edges[edgeId]->getDegree(),
            "nodeOffset of range in nodesOnEdgePinIdx");

  unsigned pOffsetIdx = _nodesOnEdgePinIdx[_noeBeginEnd[edgeId] + nodeOffset];

  abkassert2(pOffsetIdx < _pinOffsets.size(),
             "referenced pin is out of range. asked for pin: ", pOffsetIdx);

  // pinOffset is the offset if the cell has orientation 'N'
  const PinPoint& pinPoint = _pinOffsets[pOffsetIdx];
  Point pinOffset(static_cast<double>(pinPoint.offset.x),
                  static_cast<double>(pinPoint.offset.y));
  itHGFNodeLocal ndItr = _edges[edgeId]->nodesBegin() + nodeOffset;
  unsigned nodeId = (*ndItr)->getIndex();

  calcRealOffset(nodeId, pinOffset, cellOri);
  return pinOffset;
}

inline BBox HGraphWDimensions::nodesOnEdgePinOffsetBBox(unsigned nodeOffset,
                                                        unsigned edgeId,
                                                        Orient cellOri) const {
  abkassert(edgeId < _edges.size(), "edgeId out of range in nodesOnEdgePinIdx");
  abkassert(nodeOffset < _edges[edgeId]->getDegree(),
            "nodeOffset of range in nodesOnEdgePinIdx");

  unsigned pOffsetIdx = _nodesOnEdgePinIdx[_noeBeginEnd[edgeId] + nodeOffset];

  abkassert2(pOffsetIdx < _pinOffsets.size(),
             "referenced pin is out of range. asked for pin: ", pOffsetIdx);

  // pinOffset is the offset if the cell has orientation 'N'
  const PinPoint& pinPoint = _pinOffsets[pOffsetIdx];

  Point pinCenter(static_cast<double>(pinPoint.offset.x),
                  static_cast<double>(pinPoint.offset.y));
  itHGFNodeLocal ndItr = _edges[edgeId]->nodesBegin() + nodeOffset;
  unsigned nodeId = (*ndItr)->getIndex();
  calcRealOffset(nodeId, pinCenter, cellOri);

  BBox pinBBox;
  if (pinPoint.halfBBoxDims.x == 0. && pinPoint.halfBBoxDims.y == 0.) {
    pinBBox.xMin = pinBBox.xMax = pinCenter.x;
    pinBBox.yMin = pinBBox.yMax = pinCenter.y;
  } else {
    unsigned angle = cellOri.getAngle();
    if (angle == 0 || angle == 180) {
      pinBBox.xMin = pinCenter.x - static_cast<double>(pinPoint.halfBBoxDims.x);
      pinBBox.xMax = pinCenter.x + static_cast<double>(pinPoint.halfBBoxDims.x);
      pinBBox.yMin = pinCenter.y - static_cast<double>(pinPoint.halfBBoxDims.y);
      pinBBox.yMax = pinCenter.y + static_cast<double>(pinPoint.halfBBoxDims.y);
    } else {
      pinBBox.xMin = pinCenter.x - static_cast<double>(pinPoint.halfBBoxDims.y);
      pinBBox.xMax = pinCenter.x + static_cast<double>(pinPoint.halfBBoxDims.y);
      pinBBox.yMin = pinCenter.y - static_cast<double>(pinPoint.halfBBoxDims.x);
      pinBBox.yMax = pinCenter.y + static_cast<double>(pinPoint.halfBBoxDims.x);
    }
  }

  return pinBBox;
}

inline Point HGraphWDimensions::edgesOnNodePinOffset(unsigned edgeOffset,
                                                     unsigned nodeId,
                                                     Orient cellOri) const {
  abkassert(nodeId < _nodes.size(),
            "nodeId out of range in edgesOnNodePinOffset");
  abkassert(edgeOffset < _nodes[nodeId]->getDegree(),
            "edgeOffset out of range in edgesOnNodePinOffset");

  unsigned pOffsetIdx = _edgesOnNodePinIdx[_eonBeginEnd[nodeId] + edgeOffset];

  abkassert(pOffsetIdx < _pinOffsets.size(), "referenced pin is out of range");

  // pinOffset is the offset if the cell has orientation 'N'
  const PinPoint& pinPoint = _pinOffsets[pOffsetIdx];
  Point pinOffset(static_cast<double>(pinPoint.offset.x),
                  static_cast<double>(pinPoint.offset.y));

  calcRealOffset(nodeId, pinOffset, cellOri);
  return pinOffset;
}

inline BBox HGraphWDimensions::edgesOnNodePinOffsetBBox(unsigned edgeOffset,
                                                        unsigned nodeId,
                                                        Orient cellOri) const {
  abkassert(nodeId < _nodes.size(),
            "nodeId out of range in edgesOnNodePinOffset");
  abkassert(edgeOffset < _nodes[nodeId]->getDegree(),
            "edgeOffset out of range in edgesOnNodePinOffset");

  unsigned pOffsetIdx = _edgesOnNodePinIdx[_eonBeginEnd[nodeId] + edgeOffset];

  abkassert(pOffsetIdx < _pinOffsets.size(), "referenced pin is out of range");

  // pinOffset is the offset if the cell has orientation 'N'
  const PinPoint& pinPoint = _pinOffsets[pOffsetIdx];

  Point pinCenter(static_cast<double>(pinPoint.offset.x),
                  static_cast<double>(pinPoint.offset.y));
  calcRealOffset(nodeId, pinCenter, cellOri);

  BBox pinBBox;
  if (pinPoint.halfBBoxDims.x == 0. && pinPoint.halfBBoxDims.y == 0.) {
    pinBBox.xMin = pinBBox.xMax = pinCenter.x;
    pinBBox.yMin = pinBBox.yMax = pinCenter.y;
  } else {
    unsigned angle = cellOri.getAngle();
    if (angle == 0 || angle == 180) {
      pinBBox.xMin = pinCenter.x - static_cast<double>(pinPoint.halfBBoxDims.x);
      pinBBox.xMax = pinCenter.x + static_cast<double>(pinPoint.halfBBoxDims.x);
      pinBBox.yMin = pinCenter.y - static_cast<double>(pinPoint.halfBBoxDims.y);
      pinBBox.yMax = pinCenter.y + static_cast<double>(pinPoint.halfBBoxDims.y);
    } else {
      pinBBox.xMin = pinCenter.x - static_cast<double>(pinPoint.halfBBoxDims.y);
      pinBBox.xMax = pinCenter.x + static_cast<double>(pinPoint.halfBBoxDims.y);
      pinBBox.yMin = pinCenter.y - static_cast<double>(pinPoint.halfBBoxDims.x);
      pinBBox.yMax = pinCenter.y + static_cast<double>(pinPoint.halfBBoxDims.x);
    }
  }

  return pinBBox;
}

/*
  inline double HGraphWDimensions::nodesOnEdgePinLoadCap
  (unsigned nodeOffset, unsigned edgeId) const
  {
  abkassert(edgeId < _edges.size(),
  "edgeId out of range in nodesOnEdgePinIdx");
  abkassert(nodeOffset < _edges[edgeId]->getDegree(),
  "nodeOffset of range in nodesOnEdgePinIdx");

  unsigned pinIdx = _nodesOnEdgePinIdx[_noeBeginEnd[edgeId]+nodeOffset];

  abkassert2(pinIdx < _delayInfo.sinkCap.size(),
  "referenced pin is out of range. asked for pin: ", pinIdx);

  return _delayInfo.sinkCap[pinIdx];
  }
*/

#ifdef SIGNAL_DIRECTIONS
inline void HGraphWDimensions::addSrcWPin(HGFNode& node, HGFEdge& edge,
                                          PinPoint pt = PinPoint(0, 0, 0, 0)) {
  _srcsWPins.push_back(HGWDimsPin(node.getIndex(), edge.getIndex(), pt));
}

inline void HGraphWDimensions::addSnkWPin(HGFNode& node, HGFEdge& edge,
                                          PinPoint pt = PinPoint(0, 0, 0, 0)) {
  _snksWPins.push_back(HGWDimsPin(node.getIndex(), edge.getIndex(), pt));
}
#endif

inline void HGraphWDimensions::addSrcSnkWPin(HGFNode& node, HGFEdge& edge,
                                             PinPoint pt = PinPoint(0, 0, 0,
                                                                    0)) {
  _srcSnksWPins.push_back(HGWDimsPin(node.getIndex(), edge.getIndex(), pt));
}

inline Point HGraphWDimensions::getRealOffset(unsigned nodeId,
                                              const Point& pinOffset,
                                              const Orient& ori) const {
  Point offset = pinOffset;
  calcRealOffset(nodeId, offset, ori);
  return offset;
}

// hhchan. the following are added for soft block manipulations
inline void HGraphWDimensions::updateNodeSoftInfo(unsigned id, float minAR,
                                                  float maxAR) {
  bool isSoft = !equalFloat(minAR, maxAR);
  _isSoft[id] = isSoft;
  _nodeMaxAR[id] = maxAR;  // irrelevant if isSoft == false
  _nodeMinAR[id] = minAR;
  if (isSoft) {
    shapeNodeProperly(id);   // synchronize width/height
    makeNodePinsCenter(id);  // synchronize pin locations
  }
}

inline void HGraphWDimensions::shapeNodeProperly(unsigned id) {
  float currArea = _nodeWidths[id] * _nodeHeights[id];
  float geoMeanAR = sqrt(_nodeMinAR[id] * _nodeMaxAR[id]);

  _nodeWidths[id] = sqrt(currArea * geoMeanAR);
  _nodeHeights[id] = currArea / _nodeWidths[id];

  /*
     float currAR = _nodeWidths[id] / _nodeHeights[id];
     if (currAR > _nodeMaxAR[id])
     {
        _nodeWidths[id] = sqrt(currArea * _nodeMaxAR[id]);
        _nodeHeights[id] = currArea / _nodeWidths[id];
     }
     else if (currAR < _nodeMinAR[id])
     {
        _nodeWidths[id] = sqrt(currArea * _nodeMinAR[id]);
        _nodeHeights[id] = currArea / _nodeWidths[id];
     }
  */
}

inline void HGraphWDimensions::markMacrosAsSoft(float minAR, float maxAR) {
  for (unsigned int n = 0; n < getNumNodes(); n++)
    if (_isMacro[n]) updateNodeSoftInfo(n, minAR, maxAR);
}

inline void HGraphWDimensions::setupDefaultSoftNodes() {
  int numNodes = getNumNodes();
  _isSoft.resize(numNodes);
  _nodeMaxAR.resize(numNodes);
  _nodeMinAR.resize(numNodes);

  for (int n = 0; n < numNodes; n++) {
    float hardAR = _nodeWidths[n] / _nodeHeights[n];
    _nodeMaxAR[n] = hardAR;
    _nodeMinAR[n] = hardAR;
    _isSoft[n] = false;
  }
}

inline void HGraphWDimensions::makeSoftNodesPinsCenter() {
  for (unsigned int n = 0; n < getNumNodes(); n++)
    if (_isSoft[n]) makeNodePinsCenter(n);
}

inline void HGraphWDimensions::makeNodePinsCenter(unsigned id) {
  for (unsigned e = 0; e < _nodes[id]->getDegree(); e++) {
    unsigned pOffsetIdx = _edgesOnNodePinIdx[_eonBeginEnd[id] + e];
    _pinOffsets[pOffsetIdx].offset.x = _nodeWidths[id] / 2.0;
    _pinOffsets[pOffsetIdx].offset.y = _nodeHeights[id] / 2.0;
  }
}

// sadya. the following are used for sorting according to value field
// should be part of something common. Here for now
struct ctrSortDoublePair {
  double index;
  double value;
};

struct doublePair_less_mag {
  bool operator()(const ctrSortDoublePair& elem1,
                  const ctrSortDoublePair& elem2) {
    return (elem1.value < elem2.value);
  }
};

#endif
