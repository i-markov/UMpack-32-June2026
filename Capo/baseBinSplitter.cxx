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

// created by Andrew Caldwell on 10/17/99 caldwell@cs.ucla.edu

#ifdef _MSC_VER
#pragma warning(disable : 4786)
#endif

#include "baseBinSplitter.h"

#include <cfloat>
#include <cmath>
#include <cstdlib>
#include <cerrno>
#include <climits>
#include <vector>

#include "ABKCommon/abkassert.h"
#include "ABKCommon/abkfunc.h"
#include "PartEvals/netCutWBits.h"
#include "RBPlace/RBPlacement.h"

using std::cout;
using std::endl;
using std::make_pair;
using std::max;
using std::min;
using std::pair;
using std::string;
using std::vector;

static double getRealNodeWidth(unsigned index, const CapoPlacer& capo,
                               const HGraphWDimensions& hgraph);
static double getRealNodeHeight(unsigned index, const CapoPlacer& capo,
                                const HGraphWDimensions& hgraph);

namespace {

unsigned parseUnsignedString(const std::string& value, const char* context) {
  abkfatal(!value.empty(), context);
  errno = 0;
  char* end = NULL;
  const unsigned long parsed = strtoul(value.c_str(), &end, 10);
  abkfatal(end != value.c_str() && *end == '\0' && errno != ERANGE &&
               parsed <= UINT_MAX,
           context);
  return static_cast<unsigned>(parsed);
}

unsigned parseCellGroupId(const std::string& cellName, unsigned maxGroupId) {
  std::string::size_type sepPos = cellName.rfind('_');
  abkfatal(sepPos != std::string::npos, "Bad cellname");
  unsigned groupId =
      parseUnsignedString(cellName.substr(sepPos + 1),
                          "Bad cellname group suffix");
  abkfatal(groupId <= maxGroupId, "Bad groupId");
  abkfatal(groupId > 0, "Bad groupId");
  return groupId;
}

}  // namespace

// this verson finds the best row to split at
void BaseBinSplitter::createHSplitBins(bool isLarge,
                                       PartitioningProblemForCapo& problem) {
  vector<double> newBinCellAreas(2, 0.);
  vector<vector<unsigned> > cellsInNewBins(2);

  const SubHGraph& hgraph = problem.getSubHGraph();
  const Partitioning& soln = problem.getBestSoln();

  for (unsigned nId = hgraph.getNumTerminals(); nId != hgraph.getNumNodes();
       ++nId) {
    unsigned origId = hgraph.newNode2OrigIdx(nId);
    unsigned part = soln[nId].isInPart(0) ? 0 : 1;

    cellsInNewBins[part].push_back(origId);
    newBinCellAreas[part] += hgraph.getWeight(nId);
  }

  vector<double> actualCaps(2, 0.);
  double p0WSWeight = 1.;
  double p1WSWeight = 1.;
  unsigned splitRow =
      findBestSplitRow(newBinCellAreas, p0WSWeight, p1WSWeight, actualCaps);

  if (_params.verb.getForMajStats() > 4)
    cout << "Computed splitrow is " << splitRow << endl;

  // destroy the hgraph from the problem
  // and rebuild the main hgraph if necessary

  problem.destroyHGraph();

  if (isLarge) {
    // if it wasn't actually destroyed, this is a no-op
    _capo.rebuildHGraph();
    _capo.getParams().maxMem->update(
        "Split Large H after rebuilding the HGraph");
  }

  createHSplitBins(cellsInNewBins, splitRow, problem.getMaxCapacities()[0][0],
                   problem.getMaxCapacities()[1][0]);
}

// this version uses the splitRow it's given
void BaseBinSplitter::createHSplitBins(bool isLarge,
                                       PartitioningProblemForCapo& problem,
                                       unsigned splitRow) {
  vector<vector<unsigned> > cellsInNewBins(2);

  const SubHGraph& hgraph = problem.getSubHGraph();
  const Partitioning& soln = problem.getBestSoln();
  unsigned nId;
  for (nId = hgraph.getNumTerminals(); nId != hgraph.getNumNodes(); nId++) {
    unsigned origId = hgraph.newNode2OrigIdx(nId);
    unsigned part = soln[nId].isInPart(0) ? 0 : 1;

    cellsInNewBins[part].push_back(origId);
  }

  // destroy the hgraph from the problem
  // and rebuild the main hgraph if necessary

  problem.destroyHGraph();

  if (isLarge) {
    // if it wasn't actually destroyed, this is a no-op
    _capo.rebuildHGraph();
    _capo.getParams().maxMem->update(
        "Split Large H after rebuilding the HGraph");
  }

  createHSplitBins(cellsInNewBins, splitRow, problem.getMaxCapacities()[0][0],
                   problem.getMaxCapacities()[1][0]);
}

void BaseBinSplitter::createHSplitBins(
    const vector<vector<unsigned> >& cellsInNewBins, unsigned splitRow,
    double p0MaxCap, double p1MaxCap) {
  const vector<const RBPCoreRow*>& rows = _bin.getRows();
  double ySplit = rows[splitRow]->getYCoord();

  _splitRow = splitRow;

  vector<unsigned> oracleCellsP0, oracleCellsP1;
  for (vector<unsigned>::const_iterator i = _bin.getOracleCellIds().begin();
       i != _bin.getOracleCellIds().end(); ++i) {
    unsigned angle = _capo.getOraclePlacement().getOrient(*i).getAngle();
    bool notRotated = angle == 0 || angle == 180;
    double halfCellHeight =
        notRotated ? 0.5 * _capo.getNetlistHGraph().getNodeHeight(*i)
                   : 0.5 * _capo.getNetlistHGraph().getNodeWidth(*i);
    if (greaterOrEqualDouble(_capo.getOraclePlacement()[*i].y + halfCellHeight,
                             ySplit)) {
      oracleCellsP0.push_back(*i);
    } else {
      oracleCellsP1.push_back(*i);
    }
  }

  vector<double> newBinSiteAreas(2);

  for (unsigned r = 0; r < splitRow; r++)
    newBinSiteAreas[1] += _bin.getContainedAreaInCoreRow(r);

  newBinSiteAreas[0] = _bin.getCapacity() - newBinSiteAreas[1];

  //    const vector<const RBPCoreRow*>& rows = _bin.getRows();

  itCBRowOffset sBegin = _bin._startOffsets.begin();
  itCBRowOffset sEnd = _bin._startOffsets.end();
  itCBRowOffset eBegin = _bin._endOffsets.begin();
  itCBRowOffset eEnd = _bin._endOffsets.end();

  vector<CROffset> p0StartOffsets(sBegin + splitRow, sEnd);
  vector<CROffset> p0EndOffsets(eBegin + splitRow, eEnd);
  vector<CROffset> p1StartOffsets(sBegin, sBegin + splitRow);
  vector<CROffset> p1EndOffsets(eBegin, eBegin + splitRow);

  //    double ySplit = rows[splitRow]->getYCoord();

  //    _bin.unLinkNeighbors();  // will link children instead
  // determine neighbors of forthcoming child bins
  vector<CapoBin*> p0LeftNeighbors, p1LeftNeighbors;

  unsigned k;
  for (k = 0; k != _bin._leftNeighbors.size(); k++) {
    CapoBin* tempB = _bin._leftNeighbors[k];
    double yMax = tempB->getBBox().yMax;
    double yMin = tempB->getBBox().yMin;
    if (yMax > ySplit) p0LeftNeighbors.push_back(tempB);
    if (yMin < ySplit) p1LeftNeighbors.push_back(tempB);
  }

  vector<CapoBin*> p0RightNeighbors, p1RightNeighbors;
  for (k = 0; k != _bin._rightNeighbors.size(); k++) {
    CapoBin* tempB = _bin._rightNeighbors[k];
    double yMax = tempB->getBBox().yMax;
    double yMin = tempB->getBBox().yMin;
    if (yMax > ySplit) p0RightNeighbors.push_back(tempB);
    if (yMin < ySplit) p1RightNeighbors.push_back(tempB);
  }

  vector<CapoBin*> noNeighbors;  // to prevent duplicate links
                                 // between child bin... see below

  const HGraphWDimensions& nlHGraph = _capo.getNetlistHGraph();

  string binName = _bin.getName() + string("_H0");

  if (cellsInNewBins[0].size() > 0) {
    _newBins.push_back(new CapoBin(
        cellsInNewBins[0], _bin.rowsBegin() + splitRow, _bin.rowsEnd(),
        p0StartOffsets, p0EndOffsets, _bin._neighborsAbove,
        (cellsInNewBins[1].size() > 0) ? noNeighbors : _bin._neighborsBelow,
        p0LeftNeighbors, p0RightNeighbors, _capo, binName));
    _newBins.back()->linkNeighbors();
    const double& binCap = _newBins.back()->getCapacity();
    _newBins.back()->setMaxRecCellArea(min(binCap, p0MaxCap));
    _newBins.back()->setOracleCellIds(oracleCellsP0);
    _newBins.back()->setEcoAllowed(_bin.isEcoAllowed());

    // <aaronnn> prune obstacle vector for this bin,
    // depending on ySplit coord, and pass it on
    vector<unsigned> obstacleCells;
    for (vector<unsigned>::const_iterator i = _bin.getObstacleCellIds().begin();
         i != _bin.getObstacleCellIds().end(); ++i) {
      if (_capo.getPlacement()[*i].y + getRealNodeHeight(*i, _capo, nlHGraph) >=
          ySplit)
        obstacleCells.push_back(*i);
    }
    _newBins.back()->setObstacleCellIds(obstacleCells);
  }

  binName = _bin.getName() + string("_H1");

  if (cellsInNewBins[1].size() > 0) {
    _newBins.push_back(new CapoBin(
        cellsInNewBins[1], _bin.rowsBegin(), _bin.rowsBegin() + splitRow,
        p1StartOffsets, p1EndOffsets,
        (cellsInNewBins[0].size() > 0) ? noNeighbors : _bin._neighborsAbove,
        _bin._neighborsBelow, p1LeftNeighbors, p1RightNeighbors, _capo,
        binName));
    _newBins.back()->linkNeighbors();
    const double& binCap = _newBins.back()->getCapacity();
    _newBins.back()->setMaxRecCellArea(min(binCap, p1MaxCap));
    _newBins.back()->setOracleCellIds(oracleCellsP1);
    _newBins.back()->setEcoAllowed(_bin.isEcoAllowed());

    // <aaronnn> prune obstacle vector for this bin,
    // depending on ySplit coord, and pass it on
    vector<unsigned> obstacleCells;
    for (vector<unsigned>::const_iterator i = _bin.getObstacleCellIds().begin();
         i != _bin.getObstacleCellIds().end(); ++i) {
      if (_capo.getPlacement()[*i].y <= ySplit) obstacleCells.push_back(*i);
    }
    _newBins.back()->setObstacleCellIds(obstacleCells);
  }

  if (cellsInNewBins[0].size() && cellsInNewBins[1].size()) {
    _newBins[0]->_addNeighborBelow(_newBins[1]);
    _newBins[1]->_addNeighborAbove(_newBins[0]);
    _newBins[0]->_sibling = _newBins[1];
    _newBins[1]->_sibling = _newBins[0];
  }
}

void BaseBinSplitter::createVSplitBins(bool isLarge,
                                       PartitioningProblemForCapo& problem,
                                       double xSplit) {
  vector<vector<unsigned> > cellsInNewBins(2);

  const Partitioning& soln = problem.getBestSoln();
  const SubHGraph& hgraph = problem.getSubHGraph();

  for (unsigned c = hgraph.getNumTerminals(); c < hgraph.getNumNodes(); c++) {
    unsigned origId = hgraph.newNode2OrigIdx(c);
    unsigned part = soln[c].isInPart(0) ? 0 : 1;

    cellsInNewBins[part].push_back(origId);
  }

  // destroy the hgraph from the problem
  // and rebuild the main hgraph if necessary

  problem.destroyHGraph();

  if (isLarge) {
    // if it wasn't actually destroyed, this is a no-op
    _capo.rebuildHGraph();
    _capo.getParams().maxMem->update(
        "Split Large V after rebuilding the HGraph");
  }

  createVSplitBins(cellsInNewBins, xSplit, problem.getMaxCapacities()[0][0],
                   problem.getMaxCapacities()[1][0]);
}

void BaseBinSplitter::createVSplitBins(
    const vector<vector<unsigned> >& cellsInNewBins, double xSplit,
    double p0MaxCap, double p1MaxCap) {
  const vector<const RBPCoreRow*>& rows = _bin.getRows();
  _xSplit = xSplit;

  vector<unsigned> oracleCellsP0, oracleCellsP1;
  for (vector<unsigned>::const_iterator i = _bin.getOracleCellIds().begin();
       i != _bin.getOracleCellIds().end(); ++i) {
    unsigned angle = _capo.getOraclePlacement().getOrient(*i).getAngle();
    bool notRotated = angle == 0 || angle == 180;
    double halfCellWidth =
        notRotated ? 0.5 * _capo.getNetlistHGraph().getNodeWidth(*i)
                   : 0.5 * _capo.getNetlistHGraph().getNodeHeight(*i);
    if (lessOrEqualDouble(_capo.getOraclePlacement()[*i].x + halfCellWidth,
                          xSplit)) {
      oracleCellsP0.push_back(*i);
    } else {
      oracleCellsP1.push_back(*i);
    }
  }

  vector<double> siteAreas(2, 0.0);
  _bin.computePartAreas(xSplit, siteAreas);

  // construct the vectors of coreRows and [start/end]offsets

  vector<const RBPCoreRow*> p0CoreRows;
  vector<const RBPCoreRow*> p1CoreRows;
  vector<CROffset> p0StartOffsets;
  vector<CROffset> p0EndOffsets;
  vector<CROffset> p1StartOffsets;
  vector<CROffset> p1EndOffsets;
  BBox p0BBox(_bin.getBBox());
  BBox p1BBox(_bin.getBBox());
  p0BBox.xMax = xSplit;
  p1BBox.xMin = xSplit;

  for (unsigned crId = 0; crId < rows.size(); crId++) {
    const RBPCoreRow& cr = *(rows[crId]);
    const double spacing = cr.getSpacing();
    const double width = cr.getSiteWidth();

    const CROffset& sOffset = _bin._startOffsets[crId];
    const CROffset& eOffset = _bin._endOffsets[crId];
    const RBPSubRow& leftSR = cr[sOffset.first];
    const RBPSubRow& rightSR = cr[eOffset.first];

    double leftSiteEdge = leftSR.getXStart() + spacing * sOffset.second;
    /*
    double rightSiteEdge = rightSR.getXStart()+
                        spacing*(eOffset.second+1);
    */
    double rightSiteEdge =
        rightSR.getXStart() + spacing * (eOffset.second) + width;
    // eOffset.second is the last CONTAINED site, so the right
    // edge of it is the left edge of the 1st not contained site...
    // or eOffset.second+1

    if (rightSiteEdge <= xSplit) {  // core row is completely inP0
      p0CoreRows.push_back(rows[crId]);
      p0StartOffsets.push_back(sOffset);
      p0EndOffsets.push_back(eOffset);
    } else if (leftSiteEdge >= xSplit) {  // core row is completely in P1
      p1CoreRows.push_back(rows[crId]);
      p1StartOffsets.push_back(sOffset);
      p1EndOffsets.push_back(eOffset);
    } else {  // the coreRow is split
      p0CoreRows.push_back(rows[crId]);
      p1CoreRows.push_back(rows[crId]);
      p0StartOffsets.push_back(sOffset);
      p1EndOffsets.push_back(eOffset);

      unsigned srId;
      for (srId = sOffset.first; srId <= eOffset.first; srId++) {
        if (cr[srId].getXEnd() <= xSplit) continue;  // totally in P0

        if (cr[srId].getXStart() >= xSplit)  // 1st one totally in P1
        {
          p1StartOffsets.push_back(CROffset(srId, 0));

          if (srId == sOffset.first) {
            p0CoreRows.pop_back();
            p0StartOffsets.pop_back();
          } else {
            p0EndOffsets.push_back(
                CROffset(srId - 1, cr[srId - 1].getNumSites() - 1));
          }

          break;
        } else  // this subRow crosses xSplit
        {
          CROffset splitOffset(srId, UINT_MAX);
          p1StartOffsets.push_back(splitOffset);
          p0EndOffsets.push_back(splitOffset);

          double siteOffset = xSplit - cr[srId].getXStart();
          siteOffset /= spacing;

          unsigned splitSite = static_cast<unsigned>(rint(siteOffset));
          //                        unsigned splitSite0 =
          //                        static_cast<unsigned>(floor(siteOffset));
          //                        unsigned splitSite1 =
          //                        static_cast<unsigned>(ceil(siteOffset));

          if (splitSite == 0) {
            if (srId == sOffset.first) {
              p0CoreRows.pop_back();
              p0StartOffsets.pop_back();
              p0EndOffsets.pop_back();
            } else {
              p0EndOffsets.back().first = srId - 1;
              p0EndOffsets.back().second = cr[srId - 1].getNumSites() - 1;
            }
          } else
            p0EndOffsets.back().second = splitSite - 1;

          if (splitSite >= cr[srId].getNumSites()) {
            if (srId == eOffset.first) {
              p1CoreRows.pop_back();
              p1StartOffsets.pop_back();
              p1EndOffsets.pop_back();
            } else {
              p1StartOffsets.back().first = srId + 1;
              p1StartOffsets.back().second = 0;
            }
          } else
            p1StartOffsets.back().second = splitSite;

          //			if(splitSite == 0) splitSite = 1;
          //                      p1StartOffsets.back().second = splitSite;
          //                      p0EndOffsets.back().second   = splitSite-1;
          break;
        }
      }
      abkfatal(p0StartOffsets.size() == p0EndOffsets.size(),
               "p0Offset size mismatch");
      abkfatal(p1StartOffsets.size() == p1EndOffsets.size(),
               "p1Offset size mismatch");
    }
  }

  // if  cellArea < siteArea in both bins, we are fine.. quit
  // Compute the _lacking site area_ in each bin. If it's the same in the two
  // bins, quit. Otherwise, if the difference is 1 site, quit. If the difference
  // is >1, then divide by 2, to find out how many sites must be moved (and we
  // also know in which direction) Allocate the moving sites to rows evenly (if
  // we run into trouble with cell packing, we will worry about that when we see
  // how this happens)
  // By sadya Create jagged bins for better WS handling
  if ((rows.size() == 2 || rows.size() == 3) && 0) {
    const HGraphWDimensions& hgraph = _capo.getNetlistHGraph();

    const RBPCoreRow& cr = *(rows[0]);
    double siteSpacing = cr.getSpacing();
    double siteHeight = cr.getHeight();

    double p0CellArea = 0, p1CellArea = 0;
    double p0SiteArea = 0, p1SiteArea = 0;
    unsigned i;
    for (i = 0; i < cellsInNewBins[0].size(); ++i) {
      unsigned idx = cellsInNewBins[0][i];
      double cellArea =
          (ceil(hgraph.getNodeWidth(idx) / siteSpacing) * siteSpacing) *
          (ceil(hgraph.getNodeHeight(idx) / siteHeight) * siteHeight);
      p0CellArea += cellArea;
    }
    for (i = 0; i < cellsInNewBins[1].size(); ++i) {
      unsigned idx = cellsInNewBins[1][i];
      double cellArea =
          (ceil(hgraph.getNodeWidth(idx) / siteSpacing) * siteSpacing) *
          (ceil(hgraph.getNodeHeight(idx) / siteHeight) * siteHeight);
      p1CellArea += cellArea;
    }

    unsigned p0NumSites = 0, p1NumSites = 0;
    for (i = 0; i < p0CoreRows.size(); ++i) {
      CROffset startOffset = p0StartOffsets[i];
      CROffset endOffset = p0EndOffsets[i];
      const RBPCoreRow& curRow = *p0CoreRows[i];
      //        const RBPSubRow& startSR = curRow[startOffset.first];
      //	  const RBPSubRow& endSR   = curRow[endOffset.first];

      for (unsigned sr = startOffset.first; sr <= endOffset.first; sr++) {
        int sitesInThisSR = curRow[sr].getNumSites();
        if (sr == endOffset.first) sitesInThisSR = endOffset.second + 1;
        if (sr == startOffset.first) sitesInThisSR -= startOffset.second;
        p0NumSites += sitesInThisSR;
      }
    }
    for (i = 0; i < p1CoreRows.size(); ++i) {
      CROffset startOffset = p1StartOffsets[i];
      CROffset endOffset = p1EndOffsets[i];
      const RBPCoreRow& curRow = *p1CoreRows[i];
      //	  const RBPSubRow& startSR = curRow[startOffset.first];
      //	  const RBPSubRow& endSR   = curRow[endOffset.first];

      for (unsigned sr = startOffset.first; sr <= endOffset.first; sr++) {
        int sitesInThisSR = curRow[sr].getNumSites();
        if (sr == endOffset.first) sitesInThisSR = endOffset.second + 1;
        if (sr == startOffset.first) sitesInThisSR -= startOffset.second;
        p1NumSites += sitesInThisSR;
      }
    }
    p0SiteArea = p0NumSites * siteHeight * siteSpacing;
    p1SiteArea = p1NumSites * siteHeight * siteSpacing;

    double p0WS = p0SiteArea - p0CellArea;
    double p1WS = p1SiteArea - p1CellArea;
    bool reDoSplit = true;
    if (p0WS < 0 || p1WS < 0)  // special case only if -ve WS in one partition
    {
      double p0WSsites =
          ceil((fabs(p0WS) / siteHeight) / siteSpacing) * siteSpacing;
      double p1WSsites =
          ceil((fabs(p1WS) / siteHeight) / siteSpacing) * siteSpacing;
      double p0ExtraSites = 0, p1ExtraSites = 0;
      if (p0WS < 0 && p1WS < 0 && fabs(p0WSsites - p1WSsites) > 1) {
        double avg = floor((p0WSsites + p1WSsites) / 2.0);
        p0ExtraSites = avg - p0WSsites;
        p1ExtraSites = -p0ExtraSites;
      } else if (p0WS < 0 && p1WS > 0 && (p1WSsites + p0WSsites) > 1) {
        if (p0WSsites > p1WSsites) {
          p0ExtraSites = p1WSsites;
          p1ExtraSites = -p1WSsites;
        } else {
          p0ExtraSites = p0WSsites;
          p1ExtraSites = -p0WSsites;
        }
      } else if (p0WS > 0 && p1WS < 0 && (p0WSsites + p1WSsites) > 1) {
        if (p1WSsites > p0WSsites) {
          p1ExtraSites = p0WSsites;
          p0ExtraSites = -p0WSsites;
        } else {
          p1ExtraSites = p1WSsites;
          p0ExtraSites = -p1WSsites;
        }
      } else
        reDoSplit = false;

      if (reDoSplit) {
        // for now do all adds and subtracts to second row
        if (p0CoreRows.size() > 1 && p1CoreRows.size() > 1) {
          CROffset& p0StartOffset = p0StartOffsets[1];
          CROffset& p0EndOffset = p0EndOffsets[1];
          const RBPCoreRow& p0CurRow = *p0CoreRows[1];

          CROffset& p1StartOffset = p1StartOffsets[1];
          //		  CROffset& p1EndOffset = p1EndOffsets[1];
          const RBPCoreRow& p1CurRow = *p1CoreRows[1];

          if (p0CurRow.getYCoord() == p1CurRow.getYCoord() &&
              p0EndOffset.first == p1StartOffset.first)  // same subrow
          {
            const RBPSubRow& SR = p0CurRow[p0StartOffset.first];

            //		      double addedSites=0;
            if (p0ExtraSites > 0)  // add sites to P0
            {
              double curEnd =
                  SR.getXStart() + (p0StartOffset.second + 1) * siteSpacing;
              double newEnd = curEnd + p0ExtraSites * siteSpacing;
              if (newEnd >= SR.getXEnd()) {
                newEnd = SR.getXEnd() - siteSpacing;
              }
              double new0EndOffset =
                  ceil((newEnd - SR.getXStart()) / siteSpacing);
              double new1StartOffset = new0EndOffset + 1;

              p0EndOffset.second = unsigned(new0EndOffset);
              p1StartOffset.second = unsigned(new1StartOffset);
            } else  // add sites to P1
            {
              double curEnd =
                  SR.getXStart() + (p0StartOffset.second + 1) * siteSpacing;
              double newEnd = curEnd + p0ExtraSites * siteSpacing;
              if (newEnd <= SR.getXStart()) {
                newEnd = SR.getXStart() + siteSpacing;
              }
              double new0EndOffset =
                  ceil((newEnd - SR.getXStart()) / siteSpacing);
              double new1StartOffset = new0EndOffset + 1;
              p0EndOffset.second = unsigned(new0EndOffset);
              p1StartOffset.second = unsigned(new1StartOffset);
            }
          }
        }
      }
    }
  }

  //    _bin.unLinkNeighbors();  // will link children instead

  // determine neighbors of forthcoming child bins
  vector<CapoBin*> p0NeighborsAbove, p1NeighborsAbove;
  unsigned k;
  for (k = 0; k != _bin._neighborsAbove.size(); k++) {
    CapoBin* tempB = _bin._neighborsAbove[k];
    double xMax = tempB->getBBox().xMax;
    double xMin = tempB->getBBox().xMin;

    if (xMin < xSplit) p0NeighborsAbove.push_back(tempB);
    if (xMax > xSplit) p1NeighborsAbove.push_back(tempB);
  }

  vector<CapoBin*> p0NeighborsBelow, p1NeighborsBelow;
  for (k = 0; k != _bin._neighborsBelow.size(); k++) {
    CapoBin* tempB = _bin._neighborsBelow[k];
    double xMax = tempB->getBBox().xMax;
    double xMin = tempB->getBBox().xMin;

    if (xMin < xSplit) p0NeighborsBelow.push_back(tempB);
    if (xMax > xSplit) p1NeighborsBelow.push_back(tempB);
  }

  // to prevent duplicate links between child bin
  vector<CapoBin*> noNeighbors;

  _newBins.clear();

  string binName = _bin.getName() + string("_V0");

  const HGraphWDimensions& nlHGraph = _capo.getNetlistHGraph();

  if (cellsInNewBins[0].size() > 0) {
    _newBins.push_back(new CapoBin(
        cellsInNewBins[0], p0CoreRows.begin(), p0CoreRows.end(), p0StartOffsets,
        p0EndOffsets, p0NeighborsAbove, p0NeighborsBelow, _bin._leftNeighbors,
        (cellsInNewBins[1].size() > 0) ? noNeighbors : _bin._rightNeighbors,
        _capo, binName));
    _newBins.back()->linkNeighbors();
    const double& binCap = _newBins.back()->getCapacity();
    _newBins.back()->setMaxRecCellArea(min(binCap, p0MaxCap));
    _newBins.back()->setOracleCellIds(oracleCellsP0);
    _newBins.back()->setEcoAllowed(_bin.isEcoAllowed());

    // <aaronnn> prune obstacle vector for this bin,
    // depending on xSplit coord, and pass it on
    vector<unsigned> obstacleCells;
    for (vector<unsigned>::const_iterator i = _bin.getObstacleCellIds().begin();
         i != _bin.getObstacleCellIds().end(); ++i) {
      if (_capo.getPlacement()[*i].x <= xSplit) obstacleCells.push_back(*i);
    }
    _newBins.back()->setObstacleCellIds(obstacleCells);
  }

  binName = _bin.getName() + string("_V1");

  if (cellsInNewBins[1].size() > 0) {
    _newBins.push_back(new CapoBin(
        cellsInNewBins[1], p1CoreRows.begin(), p1CoreRows.end(), p1StartOffsets,
        p1EndOffsets, p1NeighborsAbove, p1NeighborsBelow,
        (cellsInNewBins[0].size() > 0) ? noNeighbors : _bin._leftNeighbors,
        _bin._rightNeighbors, _capo, binName));
    _newBins.back()->linkNeighbors();
    const double& binCap = _newBins.back()->getCapacity();
    _newBins.back()->setMaxRecCellArea(min(binCap, p1MaxCap));
    _newBins.back()->setOracleCellIds(oracleCellsP1);
    _newBins.back()->setEcoAllowed(_bin.isEcoAllowed());

    // <aaronnn> prune obstacle vector for this bin,
    // depending on xSplit coord, and pass it on
    vector<unsigned> obstacleCells;
    for (vector<unsigned>::const_iterator i = _bin.getObstacleCellIds().begin();
         i != _bin.getObstacleCellIds().end(); ++i) {
      if (_capo.getPlacement()[*i].x + getRealNodeWidth(*i, _capo, nlHGraph) >=
          xSplit)
        obstacleCells.push_back(*i);
    }
    _newBins.back()->setObstacleCellIds(obstacleCells);
  }

  if (cellsInNewBins[0].size() && cellsInNewBins[1].size()) {
    _newBins[0]->_addRightNeighbor(_newBins[1]);
    _newBins[1]->_addLeftNeighbor(_newBins[0]);
    _newBins[0]->_sibling = _newBins[1];
    _newBins[1]->_sibling = _newBins[0];
  }
}

// P0 is the top partition
// capoBin rows are sorted in INCREASING Y
// so..P0 will contain rows [_splitRow->(numrows-1)]
// and P1 rows [0->(_splitRow-1)]
//(_splitRow, then, gives the number of rows in p1..or,
// equivlently, the index of the bottom row in P0)

// WARNING: side-effect! populates actualCaps

unsigned BaseBinSplitter::findBestSplitRow(const vector<double>& cellAreas,
                                           double p0WSWeight, double p1WSWeight,
                                           vector<double>& actualCaps) const {
  if (equalDouble(cellAreas[0], 0.) || equalDouble(cellAreas[1], 0.)) {
    unsigned curRow;
    if (equalDouble(cellAreas[1], 0.)) {
      if (_params.verb.getForMajStats() > 1 || _params.verb.getForActions() > 1)
        cout << "Part 1 of H split has no cell area, giving it one third of "
                "total area"
             << endl;
      curRow = max(1u, (_bin.getNumRows() + 1) / 3);
    } else {
      if (_params.verb.getForMajStats() > 1 || _params.verb.getForActions() > 1)
        cout << "Part 0 of H split has no cell area, giving it one third of "
                "total area"
             << endl;
      curRow = min(_bin.getNumRows() - 1, (2 * _bin.getNumRows() - 1) / 3);
    }
    actualCaps[0] = 0.;
    actualCaps[1] = 0.;
    for (unsigned r = 0; r < curRow; ++r)
      actualCaps[1] += _bin.getContainedAreaInCoreRow(r);
    for (unsigned r = curRow; r < _bin.getNumRows(); ++r)
      actualCaps[0] += _bin.getContainedAreaInCoreRow(r);

    return curRow;
  }

  double totalCap = _bin.getCapacity();
  unsigned splitRow = 1;
  double p1Cap = _bin.getContainedAreaInCoreRow(0);
  double p0Cap = totalCap - p1Cap;
  double relP1WS = (p1Cap - cellAreas[1]) / p1Cap;
  double relP0WS = (p0Cap - cellAreas[0]) / p0Cap;
  double bestWSDelta = fabs(p0WSWeight * relP1WS - p1WSWeight * relP0WS);
  double bestP1Cap = p1Cap;

  for (unsigned r = 1; r < _bin.getNumRows() - 1; r++) {
    p1Cap += _bin.getContainedAreaInCoreRow(r);
    p0Cap = totalCap - p1Cap;
    relP1WS = (p1Cap - cellAreas[1]) / p1Cap;
    relP0WS = (p0Cap - cellAreas[0]) / p0Cap;
    double weightedWSDelta =
        (greaterOrEqualDouble(relP0WS, 0.) || greaterOrEqualDouble(relP1WS, 0.))
            ? fabs(p0WSWeight * relP1WS - p1WSWeight * relP0WS)
            : fabs(p0WSWeight * relP0WS - p1WSWeight * relP1WS);

    if (lessThanDouble(weightedWSDelta, bestWSDelta)) {
      bestP1Cap = p1Cap;
      bestWSDelta = weightedWSDelta;
      splitRow = r + 1;
    }
  }

  if (splitRow == 0 || splitRow > _bin.getNumRows()) {
    cout << endl << endl;
    cout << "ERROR: new bin has no rows" << endl;

    cout << "Cell Areas: " << cellAreas[0] << "  " << cellAreas[1] << endl;
    cout << "SplitRow:   " << splitRow << endl;
    cout << _bin << endl;
    abkfatal(splitRow != 0 && splitRow < _bin.getNumRows(),
             "new bin has no rows");
  }

  actualCaps = vector<double>(2, 0.0);
  actualCaps[0] = _bin.getCapacity() - bestP1Cap;
  actualCaps[1] = bestP1Cap;
  return splitRow;
}

unsigned BaseBinSplitter::shiftSplitRow(unsigned curRow,
                                        const vector<double>& cellAreas,
                                        vector<double>& actualCaps,
                                        double localWS) const {
  if (equalDouble(cellAreas[1], 0.)) {
    if (_params.verb.getForMajStats() > 1 || _params.verb.getForActions() > 1)
      cout << "Part 1 of H split has no cell area, giving it one third of "
              "total area"
           << endl;
    curRow = max(1u, (_bin.getNumRows() + 1) / 3);
  } else if (equalDouble(cellAreas[0], 0.)) {
    if (_params.verb.getForMajStats() > 1 || _params.verb.getForActions() > 1)
      cout << "Part 0 of H split has no cell area, giving it one third of "
              "total area"
           << endl;
    curRow = min(_bin.getNumRows() - 1, (2 * _bin.getNumRows() - 1) / 3);
  }

  actualCaps[0] = 0.;
  actualCaps[1] = 0.;
  for (unsigned r = 0; r < curRow; ++r)
    actualCaps[1] += _bin.getContainedAreaInCoreRow(r);
  for (unsigned r = curRow; r < _bin.getNumRows(); ++r)
    actualCaps[0] += _bin.getContainedAreaInCoreRow(r);

  if (equalDouble(cellAreas[0], 0.) || equalDouble(cellAreas[1], 0.)) {
    return curRow;
  }

  double p0WS = 1. - (cellAreas[0] / actualCaps[0]);
  double p1WS = 1. - (cellAreas[1] / actualCaps[1]);

  if (greaterOrEqualDouble(p0WS, localWS) &&
      greaterOrEqualDouble(p1WS, localWS)) {
    return curRow;
  }

  double rowArea;
  if (lessThanDouble(p1WS,
                     p0WS))  // increment the split row == add more sites to p1
  {
    for (; curRow < _bin.getNumRows(); ++curRow) {
      rowArea = _bin.getContainedAreaInCoreRow(curRow);

      if (lessOrEqualDouble(
              actualCaps[0],
              rowArea)) {  // no legal solution possible so just stop now
        if (_params.verb.getForActions() > 1)
          abkwarn(0, "Forced to return an illegal row split wrt localWS");
        return curRow;
      }

      actualCaps[0] -= rowArea;
      actualCaps[1] += rowArea;

      p0WS = 1. - (cellAreas[0] / actualCaps[0]);
      p1WS = 1. - (cellAreas[1] / actualCaps[1]);

      if (greaterOrEqualDouble(p0WS, localWS) &&
          greaterOrEqualDouble(p1WS,
                               localWS)) {  // first legal solution, take it
        return (curRow + 1);
      } else if (lessThanDouble(p0WS, localWS) &&
                 greaterOrEqualDouble(
                     p1WS,
                     localWS)) {  // no legal solution possible so just stop now
        if (_params.verb.getForActions() > 1)
          abkwarn(0, "Forced to return an illegal row split wrt localWS");
        actualCaps[0] += rowArea;
        actualCaps[1] -= rowArea;
        return curRow;
      }
    }
  } else  // decrement the split row == add more sites to p0
  {
    for (; curRow > 0; --curRow) {
      rowArea = _bin.getContainedAreaInCoreRow(curRow - 1);

      if (lessOrEqualDouble(
              actualCaps[1],
              rowArea)) {  // no legal solution possible so just stop now
        if (_params.verb.getForActions() > 1)
          abkwarn(0, "Forced to return an illegal row split wrt localWS");
        return curRow;
      }

      actualCaps[0] += rowArea;
      actualCaps[1] -= rowArea;

      p0WS = 1. - (cellAreas[0] / actualCaps[0]);
      p1WS = 1. - (cellAreas[1] / actualCaps[1]);

      if (greaterOrEqualDouble(p0WS, localWS) &&
          greaterOrEqualDouble(p1WS,
                               localWS)) {  // first legal solution, take it
        return (curRow - 1);
      } else if (greaterOrEqualDouble(p0WS, localWS) &&
                 lessThanDouble(
                     p1WS,
                     localWS)) {  // no legal solution possible so just stop now
        if (_params.verb.getForActions() > 1)
          abkwarn(0, "Forced to return an illegal row split wrt localWS");
        actualCaps[0] -= rowArea;
        actualCaps[1] += rowArea;
        return curRow;
      }
    }
  }
  abkfatal(0, "Couldn't return a row split");
  return curRow;  // [MAR] to make the compiler happy
}

// splitDir=0 is horiz split    splitDir=1 is vertical split
bool BaseBinSplitter::partAreasForFillerCells(bool splitDir,
                                              vector<double>& partAreas) {
  double stdCellArea = 0;
  double fillerCellArea = 0;
  const BBox& binBBox = _bin.getBBox();
  const RBPlacement& rbplace = _capo.getRBP();
  const BBox& coreBBox = const_cast<RBPlacement&>(rbplace).getBBox(false);
  const HGraphWDimensions& hgraph = _capo.getRBP().getHGraph();
  double coreXMid = (coreBBox.xMax + coreBBox.xMin) / 2.0;
  double coreYMid = (coreBBox.yMax + coreBBox.yMin) / 2.0;
  double binXMid = (binBBox.xMax + binBBox.xMin) / 2.0;
  double binYMid = (binBBox.yMax + binBBox.yMin) / 2.0;

  const vector<unsigned>& cellIds = _bin.getCellIds();
  for (unsigned i = 0; i < cellIds.size(); ++i) {
    unsigned idx = cellIds[i];
    HGFNode& node = const_cast<HGFNode&>(hgraph.getNodeByIdx(idx));
    double cellArea = hgraph.getNodeWidth(idx) * hgraph.getNodeHeight(idx);
    unsigned numEdges = node.getNumEdges();
    if (numEdges > 0)
      stdCellArea += cellArea;
    else
      fillerCellArea += cellArea;
  }
  double totalArea = stdCellArea + fillerCellArea;
  double fillerCellRatio = fillerCellArea / totalArea;
  // double stdCellRatio = 1 - fillerCellRatio;
  if (fillerCellRatio < 0.1 || fillerCellRatio > 0.95) return false;
  double fillerCellPartRatio = fillerCellArea / (5.0 * totalArea);
  double totalCellArea = _bin.getTotalCellArea();

  bool fillerCellPart = 0;  // in which partition do filler cells go
  if (splitDir == 0)        // horiz split
  {
    if (binBBox.yMax <= coreYMid)
      fillerCellPart = 1;
    else if (binBBox.yMin >= coreYMid)
      fillerCellPart = 0;
    else if (binYMid <= coreYMid)
      fillerCellPart = 1;
    else
      fillerCellPart = 0;
  } else  // vertical split
  {
    if (binBBox.xMax <= coreXMid)
      fillerCellPart = 0;
    else if (binBBox.xMin >= coreXMid)
      fillerCellPart = 1;
    else if (binXMid <= coreXMid)
      fillerCellPart = 0;
    else
      fillerCellPart = 1;
  }

  if (fillerCellPart == 0) {
    partAreas[0] = fillerCellPartRatio * totalCellArea;
    partAreas[1] = totalCellArea - partAreas[0];
  } else {
    partAreas[1] = fillerCellPartRatio * totalCellArea;
    partAreas[0] = totalCellArea - partAreas[1];
  }
  return true;
}

static double computeWSShrinkFactorByBinDensityRatio(double densityRatio)
// <aaronnn> squeeze WS to [0.8%, 0.2%], depending on perimeter ratio from [3x,
// inf]
{
  double shrinkFactor = (-0.003 * densityRatio * densityRatio) + 0.8;
  shrinkFactor = max(shrinkFactor, 0.2);

  // double shrinkFactor = (-0.002 * densityRatio * densityRatio) + 0.8;
  // shrinkFactor = max(shrinkFactor, 0.4);

  return shrinkFactor;
}

double BaseBinSplitter::adjustSplitByUtilizationDensity(
    const PartitioningProblemForCapo& problem, double splitPt, bool isHCut)
// <aaronnn> check the partitioning solution and shift cutline if necessary, to
// reallocate bin capacities
{
  vector<double> newBinCellAreas(2, 0.0);
  vector<double> newBinCellPerims(2, 0.0);
  vector<vector<unsigned> > cellsInNewBins(2);
  vector<vector<unsigned> > obstaclesInNewBins(2);
  vector<double> binCapacities(2, 0.);
  vector<double> largestCellDimension(2, 0.);

  BBox bin0BBox;
  if (isHCut) {
    bin0BBox.add(_bin.getBBox().xMin, splitPt);
    bin0BBox.add(_bin.getBBox().xMax, _bin.getBBox().yMax);
  } else {
    bin0BBox.add(_bin.getBBox().xMin, _bin.getBBox().yMin);
    bin0BBox.add(splitPt, _bin.getBBox().yMax);
  }

  BBox bin1BBox;
  if (isHCut) {
    bin1BBox.add(_bin.getBBox().xMin, _bin.getBBox().yMin);
    bin1BBox.add(_bin.getBBox().xMax, splitPt);
  } else {
    bin1BBox.add(splitPt, _bin.getBBox().yMin);
    bin1BBox.add(_bin.getBBox().xMax, _bin.getBBox().yMax);
  }

  binCapacities[0] = _bin.getContainedAreaInBBox(bin0BBox);
  binCapacities[1] = _bin.getCapacity() - binCapacities[0];

  double minSiteSpacing = DBL_MAX;
  double minRowHeight = DBL_MAX;
  for (unsigned r = 0; r < _bin.getRows().size(); r++) {
    const RBPCoreRow& cr = *_bin.getRows()[r];
    minSiteSpacing = min(minSiteSpacing, cr.getSpacing());
    minRowHeight = min(minRowHeight, cr.getHeight());
  }

  // find cells in each bin
  const SubHGraph& hgraph = problem.getSubHGraph();
  const HGraphWDimensions& origHGraph = _capo.getRBP().getHGraph();
  const Partitioning& soln = problem.getBestSoln();
  unsigned nId;
  for (nId = hgraph.getNumTerminals(); nId != hgraph.getNumNodes(); nId++) {
    unsigned origId = hgraph.newNode2OrigIdx(nId);
    unsigned part = soln[nId].isInPart(0) ? 0 : 1;

    cellsInNewBins[part].push_back(origId);
    newBinCellAreas[part] += hgraph.getWeight(nId);
    newBinCellPerims[part] += (2 * origHGraph.getNodeWidth(origId)) +
                              (2 * origHGraph.getNodeHeight(origId));

    largestCellDimension[part] =
        max(largestCellDimension[part],
            static_cast<double>(origHGraph.getNodeWidth(origId)));
    largestCellDimension[part] =
        max(largestCellDimension[part],
            static_cast<double>(origHGraph.getNodeHeight(origId)));
  }

  // const vector<const RBPCoreRow*>& rows = _bin.getRows();

  // find obstacles in each bin
  for (vector<unsigned>::const_iterator i = _bin.getObstacleCellIds().begin();
       i != _bin.getObstacleCellIds().end(); ++i) {
    Point obstacleLoc(_capo.getPlacement()[*i].x, _capo.getPlacement()[*i].y);

    double obstacleWidth = getRealNodeWidth(*i, _capo, origHGraph);
    double obstacleHeight = getRealNodeHeight(*i, _capo, origHGraph);

    BBox obstacleBBox;
    obstacleBBox.add(obstacleLoc.x, obstacleLoc.y);
    obstacleBBox.add(obstacleLoc.x + obstacleWidth,
                     obstacleLoc.y + obstacleHeight);

    // if (obstacleLoc.y <= ySplit && obstacleLoc.y + obstacleHeight >= ySplit)
    //{
    //	// obstacle is in both bins, don't try to adjust cut to accomodate it
    //	continue
    // }
    // BBox includedObstacleBBox = obstacleBBox.intersect(_bin.getBBox());

    // unsigned part = (includedObstacleBBox.yMax >= ySplit) ? 0 : 1;

    // obstaclesInNewBins[part].push_back(*i);
    // newBinCellAreas[part] += includedObstacleBBox.getWidth()
    //     * includedObstacleBBox.getHeight();
    // newBinCellPerims[part] += (2 * includedObstacleBBox.getWidth())
    //     + (2 * includedObstacleBBox.getHeight());

    BBox includedObstacleBBox0 = obstacleBBox.intersect(bin0BBox);
    BBox includedObstacleBBox1 = obstacleBBox.intersect(bin1BBox);

    if (!includedObstacleBBox0.isEmpty()) {
      obstaclesInNewBins[0].push_back(*i);
      double obstacleArea =
          includedObstacleBBox0.getWidth() * includedObstacleBBox0.getHeight();
      newBinCellAreas[0] += obstacleArea;
      newBinCellPerims[0] += (2 * includedObstacleBBox0.getWidth()) +
                             (2 * includedObstacleBBox0.getHeight());
    }
    if (!includedObstacleBBox1.isEmpty()) {
      obstaclesInNewBins[1].push_back(*i);
      double obstacleArea =
          includedObstacleBBox1.getWidth() * includedObstacleBBox1.getHeight();
      newBinCellAreas[1] += obstacleArea;
      newBinCellPerims[1] += (2 * includedObstacleBBox1.getWidth()) +
                             (2 * includedObstacleBBox1.getHeight());
    }
  }

  // calculate the desired WS

  double bin0WS = binCapacities[0] - newBinCellAreas[0];
  double bin1WS = binCapacities[1] - newBinCellAreas[1];
  // double bin0WSPct = bin0WS / binCapacities[0];
  // double bin1WSPct = bin1WS / binCapacities[1];

  // if (cellsInNewBins[0].size() == 1 && obstaclesInNewBins[0].size() == 0
  //     && cellsInNewBins[1].size() > 20)
  //{
  //	// give just enough space for the single cell
  //	unsigned cellId = cellsInNewBins[0][0];
  //	double nodeWidth = origHGraph.getNodeWidth(cellId);
  //	double nodeHeight = origHGraph.getNodeHeight(cellId);

  //	double largestDimension = max(nodeWidth, nodeHeight);

  //	if (isHCut)
  //		return bin0BBox.yMax - min(largestDimension,
  //bin0BBox.getHeight()); 	else 		return bin0BBox.xMin + min(largestDimension,
  //bin0BBox.getWidth());
  //}
  // else if (cellsInNewBins[1].size() == 1 && obstaclesInNewBins[1].size() == 0
  //    && cellsInNewBins[0].size() > 20)
  //{
  //	// give just enough space for the single cell
  //	unsigned cellId = cellsInNewBins[1][0];
  //	double nodeWidth = origHGraph.getNodeWidth(cellId);
  //	double nodeHeight = origHGraph.getNodeHeight(cellId);

  //	double largestDimension = max(nodeWidth, nodeHeight);

  //	if (isHCut)
  //		return bin0BBox.yMin + min(largestDimension,
  //bin1BBox.getHeight()); 	else 		return bin0BBox.xMax - min(largestDimension,
  //bin1BBox.getWidth());
  //}
  if (newBinCellPerims[1] > 2 * newBinCellPerims[0] &&
      newBinCellPerims[0] > 0.0001 &&
      (obstaclesInNewBins[0].size() ==
       0 /*|| cellsInNewBins[0].size() == 1*/)) {
    // partition 0 is denser than partition 1 - give more WS to partition 1
    if (newBinCellAreas[0] > binCapacities[0] ||
        newBinCellAreas[1] > binCapacities[1]) {
      // cout << "ANDBG returning splitPt unmodified " << endl;
      return splitPt;
    }

    // double newBin0WSPct = bin0WSPct * 0.5;
    // double newBinArea = newBinCellAreas[0] + (newBin0WSPct *
    // binCapacities[0]); double newBinArea = newBinCellAreas[0] + (bin0WS *
    // 0.5);

    // squeeze WS to [0.8%, 0.2%], depending on perimeter ratio from [3x, inf]
    double newBinPerimRatio = newBinCellPerims[1] / newBinCellPerims[0];
    double newBin0WSPctFactor =
        computeWSShrinkFactorByBinDensityRatio(newBinPerimRatio);
    // cout << "**************** ANDBG 0 ratio: " << newBinCellPerims[1] /
    // newBinCellPerims[0] << endl; cout << "**************** ANDBG 0
    // wspctfactor: " << newBin0WSPctFactor << endl;
    double newBin0WSPct = newBin0WSPctFactor * bin0WS;
    double newBinArea = newBinCellAreas[0] + newBin0WSPct;

    if (isHCut) {
      double newBinWidth = bin0BBox.getWidth();
      double newBinHeight = newBinArea / newBinWidth;
      double newSplitPt = bin0BBox.yMax - newBinHeight;

      // if (newBinHeight < largestCellDimension[0])
      //	newSplitPt = max(bin0BBox.yMax - largestCellDimension[0],
      //splitPt);

      if (fabs(newSplitPt - splitPt) < minRowHeight * 2. ||
          fabs(newSplitPt - splitPt) < 0.05 * _bin.getBBox().getHeight())
        return splitPt;

      // cout << "***************************************** ANDBG shifting horiz
      // cut up from " << splitPt << " to "
      //	<< newSplitPt << endl;

      return newSplitPt;
    } else {
      double newBinHeight = bin0BBox.getHeight();
      double newBinWidth = newBinArea / newBinHeight;
      double newSplitPt = bin0BBox.xMin + newBinWidth;

      // if (newBinWidth < largestCellDimension[0])
      //	newSplitPt = min(bin0BBox.xMin + largestCellDimension[0],
      //splitPt);

      if (fabs(newSplitPt - splitPt) < minSiteSpacing * 4. ||
          fabs(newSplitPt - splitPt) < 0.05 * _bin.getBBox().getWidth())
        return splitPt;

      // cout << "***************************************** ANDBG shifting vert
      // cut left from " << splitPt << " to "
      //	<< newSplitPt << endl;

      return newSplitPt;
    }
  } else if (newBinCellPerims[0] > 2 * newBinCellPerims[1] &&
             newBinCellPerims[1] > 0.0001 &&
             (obstaclesInNewBins[1].size() ==
              0 /*|| cellsInNewBins[1].size() == 1*/)) {
    // partition 1 is denser than partition 0 - give more WS to partition 0
    if (newBinCellAreas[0] > binCapacities[0] ||
        newBinCellAreas[1] > binCapacities[1]) {
      // cout << "ANDBG returning splitPt unmodified " << endl;
      return splitPt;
    }

    // double newBin1WSPct = bin1WSPct * 0.5;
    // double newBinArea = newBinCellAreas[1] + (newBin1WSPct *
    // binCapacities[1]); double newBinArea = newBinCellAreas[1] + (bin1WS *
    // 0.5);

    // squeeze WS to [0.8%, 0.2%], depending on perimeter ratio from [3x, inf]
    double newBinPerimRatio = newBinCellPerims[0] / newBinCellPerims[1];
    double newBin1WSPctFactor =
        computeWSShrinkFactorByBinDensityRatio(newBinPerimRatio);
    // cout << "**************** ANDBG 1 ratio: " << newBinCellPerims[0] /
    // newBinCellPerims[1] << endl; cout << "**************** ANDBG 1
    // wspctfactor: " << newBin1WSPctFactor << endl;
    double newBin1WSPct = newBin1WSPctFactor * bin1WS;
    double newBinArea = newBinCellAreas[1] + newBin1WSPct;

    if (isHCut) {
      double newBinWidth = bin1BBox.getWidth();
      double newBinHeight = newBinArea / newBinWidth;
      double newSplitPt = bin1BBox.yMin + newBinHeight;

      // if (newBinHeight < largestCellDimension[1])
      //	newSplitPt = min(bin0BBox.yMin + largestCellDimension[0],
      //splitPt);

      if (fabs(newSplitPt - splitPt) < minRowHeight * 2. ||
          fabs(newSplitPt - splitPt) < 0.05 * _bin.getBBox().getHeight())
        return splitPt;

      // cout << "***************************************** ANDBG shifting horiz
      // cut down from " << splitPt << " to "
      //	<< newSplitPt << endl;

      return newSplitPt;
    } else {
      double newBinHeight = bin1BBox.getHeight();
      double newBinWidth = newBinArea / newBinHeight;
      double newSplitPt = bin1BBox.xMax - newBinWidth;

      // if (newBinWidth < largestCellDimension[1])
      //	newSplitPt = max(bin0BBox.xMax - largestCellDimension[0],
      //splitPt);

      if (fabs(newSplitPt - splitPt) < minSiteSpacing * 4. ||
          fabs(newSplitPt - splitPt) < 0.05 * _bin.getBBox().getWidth())
        return splitPt;

      // cout << "***************************************** ANDBG shifting vert
      // cut right from " << splitPt << " to "
      //	<< newSplitPt << endl;

      return newSplitPt;
    }
  }

  // cout << "ANDBG returning splitPt unmodified " << endl;
  return splitPt;
}

double getRealNodeWidth(unsigned index, const CapoPlacer& capo,
                        const HGraphWDimensions& hgraph) {
  unsigned angle = capo.getPlacement().getOrient(index).getAngle();
  bool notRotated = angle == 0 || angle == 180;
  double nodeWidth =
      notRotated ? hgraph.getNodeWidth(index) : hgraph.getNodeHeight(index);

  return nodeWidth;
}

double getRealNodeHeight(unsigned index, const CapoPlacer& capo,
                         const HGraphWDimensions& hgraph) {
  unsigned angle = capo.getPlacement().getOrient(index).getAngle();
  bool notRotated = angle == 0 || angle == 180;
  double nodeHeight =
      notRotated ? hgraph.getNodeHeight(index) : hgraph.getNodeWidth(index);

  return nodeHeight;
}

pair<double, double> BaseBinSplitter::getCongestionRatiosH(
    unsigned splitRow, const ISPD04CongMap& congmap) const {
  // p0 is the top
  BBox p0BBox = _bin.getBBox();
  p0BBox.yMin = _bin.getRows()[splitRow]->getYCoord();
  pair<double, double> p0Demand = congmap.getDemand(p0BBox);
  pair<double, double> p0Supply = congmap.getSupply(p0BBox);

  // p1 is the bottom
  BBox p1BBox = _bin.getBBox();
  p1BBox.yMax = _bin.getRows()[splitRow]->getYCoord();
  pair<double, double> p1Demand = congmap.getDemand(p1BBox);
  pair<double, double> p1Supply = congmap.getSupply(p1BBox);

  pair<double, double> ratio;

  //  both demand
  //  ratio =
  //  make_pair(p0Demand.first+p0Demand.second,p1Demand.first+p1Demand.second);

  //  hv demand
  //  ratio = make_pair(p0Demand.first,p1Demand.first);

  //  both demand ratio
  ratio = make_pair(
      (p0Demand.first + p0Demand.second) / (p0Supply.first + p0Supply.second),
      (p1Demand.first + p1Demand.second) / (p1Supply.first + p1Supply.second));

  //  hv demand ratio
  //  ratio =
  //  make_pair(p0Demand.first/p0Supply.first,p1Demand.first/p1Supply.first);

  //  vh demand ratio
  //  ratio =
  //  make_pair(p0Demand.second/p0Supply.second,p1Demand.second/p1Supply.second);

  double minpart = min(ratio.first, ratio.second);
  if (greaterThanDouble(minpart, 0.)) {
    ratio.first = ratio.first / minpart;
    ratio.second = ratio.second / minpart;
  }

  ratio.first = pow(ratio.first, _params.congExponent);
  ratio.second = pow(ratio.second, _params.congExponent);

  return ratio;
}

bool BaseBinSplitter::shouldCongestionShiftH(unsigned splitRow,
                                             const ISPD04CongMap& congmap,
                                             double avgCong) const {
  if (_params.alwaysCongShift) return true;

  // p0 is the top
  BBox p0BBox = _bin.getBBox();
  p0BBox.yMin = _bin.getRows()[splitRow]->getYCoord();
  pair<double, double> p0Demand = congmap.getDemand(p0BBox);

  // p1 is the bottom
  BBox p1BBox = _bin.getBBox();
  p1BBox.yMax = _bin.getRows()[splitRow]->getYCoord();
  pair<double, double> p1Demand = congmap.getDemand(p1BBox);

  double p0AvgCong = (p0Demand.first + p0Demand.second) / p0BBox.getArea();
  double p1AvgCong = (p1Demand.first + p1Demand.second) / p1BBox.getArea();

  return (greaterThanDouble(p0AvgCong, avgCong) ||
          greaterThanDouble(p1AvgCong, avgCong));
}

pair<double, double> BaseBinSplitter::getCongestionRatiosV(
    double xSplit, const ISPD04CongMap& congmap) const {
  // p0 is the left
  BBox p0BBox = _bin.getBBox();
  p0BBox.xMax = xSplit;
  pair<double, double> p0Demand = congmap.getDemand(p0BBox);
  pair<double, double> p0Supply = congmap.getSupply(p0BBox);

  // p1 is the right
  BBox p1BBox = _bin.getBBox();
  p1BBox.xMin = xSplit;
  pair<double, double> p1Demand = congmap.getDemand(p1BBox);
  pair<double, double> p1Supply = congmap.getSupply(p1BBox);

  pair<double, double> ratio;

  //  both demand
  //  ratio =
  //  make_pair(p0Demand.first+p0Demand.second,p1Demand.first+p1Demand.second);

  //  hv demand
  //  ratio = make_pair(p0Demand.second,p1Demand.second);

  //  both demand ratio
  ratio = make_pair(
      (p0Demand.first + p0Demand.second) / (p0Supply.first + p0Supply.second),
      (p1Demand.first + p1Demand.second) / (p1Supply.first + p1Supply.second));

  //  hv demand ratio
  //  ratio =
  //  make_pair(p0Demand.second/p0Supply.second,p1Demand.second/p1Supply.second);

  //  vh demand ratio
  //  ratio =
  //  make_pair(p0Demand.first/p0Supply.first,p1Demand.first/p1Supply.first);

  double minpart = min(ratio.first, ratio.second);
  if (greaterThanDouble(minpart, 0.)) {
    ratio.first = ratio.first / minpart;
    ratio.second = ratio.second / minpart;
  }

  ratio.first = pow(ratio.first, _params.congExponent);
  ratio.second = pow(ratio.second, _params.congExponent);

  return ratio;
}

bool BaseBinSplitter::shouldCongestionShiftV(double xSplit,
                                             const ISPD04CongMap& congmap,
                                             double avgCong) const {
  if (_params.alwaysCongShift) return true;

  // p0 is the left
  BBox p0BBox = _bin.getBBox();
  p0BBox.xMax = xSplit;
  pair<double, double> p0Demand = congmap.getDemand(p0BBox);

  // p1 is the right
  BBox p1BBox = _bin.getBBox();
  p1BBox.xMin = xSplit;
  pair<double, double> p1Demand = congmap.getDemand(p1BBox);

  double p0AvgCong = (p0Demand.first + p0Demand.second) / p0BBox.getArea();
  double p1AvgCong = (p1Demand.first + p1Demand.second) / p1BBox.getArea();

  return (greaterThanDouble(p0AvgCong, avgCong) ||
          greaterThanDouble(p1AvgCong, avgCong));
}

void BaseBinSplitter::shiftOraclePlacementH(double ecoSplitRowCoord,
                                            double shiftedSplitRowCoord) {
  double oldCutlineDistanceBot = ecoSplitRowCoord - _bin.getBBox().yMin;
  double newCutlineDistanceBot = shiftedSplitRowCoord - _bin.getBBox().yMin;
  double oldCutlineDistanceTop = _bin.getBBox().yMax - ecoSplitRowCoord;
  double newCutlineDistanceTop = _bin.getBBox().yMax - shiftedSplitRowCoord;
  double ratioBot = newCutlineDistanceBot / oldCutlineDistanceBot;
  double ratioTop = newCutlineDistanceTop / oldCutlineDistanceTop;

  for (unsigned i = 0; i < _bin.getOracleCellIds().size(); ++i) {
    unsigned cellId = _bin.getOracleCellIds()[i];
    unsigned angle = _capo.getOraclePlacement().getOrient(cellId).getAngle();
    bool notRotated = angle == 0 || angle == 180;
    double halfCellHeight =
        notRotated ? 0.5 * _capo.getNetlistHGraph().getNodeHeight(cellId)
                   : 0.5 * _capo.getNetlistHGraph().getNodeWidth(cellId);
    double cellYCenter = _capo.getOraclePlacement()[cellId].y + halfCellHeight;

    if (greaterOrEqualDouble(cellYCenter, ecoSplitRowCoord)) {
      double oldYdistTop = _bin.getBBox().yMax - cellYCenter;
      double newYcoord =
          _bin.getBBox().yMax - oldYdistTop * ratioTop - halfCellHeight;
      _capo.getMutableOraclePlacement()[cellId].y = newYcoord;
    } else {
      double oldYdistBot = cellYCenter - _bin.getBBox().yMin;
      double newYcoord =
          _bin.getBBox().yMin + oldYdistBot * ratioBot - halfCellHeight;
      _capo.getMutableOraclePlacement()[cellId].y = newYcoord;
    }
  }
}

void BaseBinSplitter::shiftOraclePlacementV(double ecoXSplit,
                                            double shiftedXSplit) {
  double oldCutlineDistanceLeft = ecoXSplit - _bin.getBBox().xMin;
  double newCutlineDistanceLeft = shiftedXSplit - _bin.getBBox().xMin;
  double oldCutlineDistanceRight = _bin.getBBox().xMax - ecoXSplit;
  double newCutlineDistanceRight = _bin.getBBox().xMax - shiftedXSplit;
  double ratioLeft = newCutlineDistanceLeft / oldCutlineDistanceLeft;
  double ratioRight = newCutlineDistanceRight / oldCutlineDistanceRight;

  for (unsigned i = 0; i < _bin.getOracleCellIds().size(); ++i) {
    unsigned cellId = _bin.getOracleCellIds()[i];
    unsigned angle = _capo.getOraclePlacement().getOrient(cellId).getAngle();
    bool notRotated = angle == 0 || angle == 180;
    double halfCellWidth =
        notRotated ? 0.5 * _capo.getNetlistHGraph().getNodeWidth(cellId)
                   : 0.5 * _capo.getNetlistHGraph().getNodeHeight(cellId);
    double cellXCenter = _capo.getOraclePlacement()[cellId].x + halfCellWidth;

    if (lessOrEqualDouble(cellXCenter, ecoXSplit)) {
      double oldXdistLeft = cellXCenter - _bin.getBBox().xMin;
      double newXcoord =
          _bin.getBBox().xMin + oldXdistLeft * ratioLeft - halfCellWidth;
      _capo.getMutableOraclePlacement()[cellId].x = newXcoord;
    } else {
      double oldXdistRight = _bin.getBBox().xMax - cellXCenter;
      double newXcoord =
          _bin.getBBox().xMax - oldXdistRight * ratioRight - halfCellWidth;
      _capo.getMutableOraclePlacement()[cellId].x = newXcoord;
    }
  }
}

bool BaseBinSplitter::shouldPartitionByGroups(
    const PartitioningProblemForCapo& problem) const {
  if (_bin.getNumCells() > _params.cellGroupsCutoff) return false;
  // check to see if we should partition this problem using groups
  unsigned numGroups = 0;
  if (_params.numCellGroups > 1) {
    // check to see how many groups there are among the cells we have to
    // partition
    BitBoard actualGroups(_params.numCellGroups);
    for (unsigned n = problem.getHGraph().getNumTerminals();
         n < problem.getHGraph().getNumNodes(); ++n) {
      unsigned origId = problem.getSubHGraph().newNode2OrigIdx(n);
      std::string cellName =
          _capo.getNetlistHGraph().getNodeNameByIndex(origId).c_str();
      unsigned groupId = parseCellGroupId(cellName, _params.numCellGroups);
      actualGroups.setBit(groupId);
    }
    numGroups = actualGroups.getIndicesOfSetBits().size();
  }

  if (numGroups <= 1) {
    // don't partition by groups, return false
    return false;
  } else {
    return true;
  }
}

double BaseBinSplitter::partitionByGroups(PartitioningProblemForCapo& problem) {
  if (_bin.getNumCells() > _params.cellGroupsCutoff) return -1.;
  // check to see if we should partition this problem using groups
  std::vector<unsigned> groupNums;
  std::vector<unsigned> movableGroupIds(problem.getHGraph().getNumNodes(),
                                        UINT_MAX);
  if (_params.numCellGroups > 1) {
    // check to see how many groups there are among the cells we have to
    // partition
    BitBoard actualGroups(_params.numCellGroups);
    for (unsigned n = problem.getHGraph().getNumTerminals();
         n < problem.getHGraph().getNumNodes(); ++n) {
      unsigned origId = problem.getSubHGraph().newNode2OrigIdx(n);
      std::string cellName =
          _capo.getNetlistHGraph().getNodeNameByIndex(origId).c_str();
      unsigned groupId = parseCellGroupId(cellName, _params.numCellGroups);
      movableGroupIds[n] = groupId;
      actualGroups.setBit(groupId);
    }
    groupNums = actualGroups.getIndicesOfSetBits();
  }

  if (groupNums.size() <= 1) {
    // don't partition by groups, return sentinel
    return -1.;
  } else {
    // do partition by groups, return true when done
    Partitioning soln(problem.getHGraph().getNumNodes());
    Partitioning best(problem.getHGraph().getNumNodes());

    for (unsigned n = 0; n < problem.getHGraph().getNumTerminals(); ++n) {
      unsigned fixed = problem.getFixedConstr()[n].getUnsigned();
      soln[n].loadBitsFrom(fixed & 3u);
    }

    double bestCut = DBL_MAX;

    for (unsigned i = 0; i < groupNums.size(); ++i) {
      // pull out one of the groups and put it in part 0
      for (unsigned n = problem.getHGraph().getNumTerminals();
           n < problem.getHGraph().getNumNodes(); ++n) {
        unsigned groupId = movableGroupIds[n];
        if (groupId == groupNums[i]) {
          soln[n].setToPart(0);
        } else {
          soln[n].setToPart(1);
        }
      }
      // evaluate the cut of this partitionment
      NetCutWBits cutEval(problem, soln);
      double cut = 0.;
      for (unsigned e = 0; e < problem.getHGraph().getNumEdges(); ++e) {
        cut += cutEval.getNetCostDouble(e) *
               problem.getHGraph().getEdgeByIdx(e).getWeight();
      }
      if (cut < bestCut) {
        bestCut = cut;
        best = soln;
      }

      if (groupNums.size() > 2) {
        // pull out one of the groups and put it in part 1
        for (unsigned n = problem.getHGraph().getNumTerminals();
             n < problem.getHGraph().getNumNodes(); ++n) {
          unsigned groupId = movableGroupIds[n];
          if (groupId == groupNums[i]) {
            soln[n].setToPart(1);
          } else {
            soln[n].setToPart(0);
          }
        }
        // evaluate the cut of this partitionment
        NetCutWBits cutEval2(problem, soln);
        cut = 0.;
        for (unsigned e = 0; e < problem.getHGraph().getNumEdges(); ++e) {
          cut += cutEval2.getNetCostDouble(e) *
                 problem.getHGraph().getEdgeByIdx(e).getWeight();
        }
        if (cut < bestCut) {
          bestCut = cut;
          best = soln;
        }
      }
    }

    problem.reserveBuffers(1);
    problem.getSolnBuffers()[0] = best;
    problem.setBestSolnNum(0);

    return bestCut;
  }
}
