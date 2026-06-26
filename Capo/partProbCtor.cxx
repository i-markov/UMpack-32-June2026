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

#ifdef _MSC_VER
#pragma warning(disable : 4786)
#endif

#include "Partitioning/termiProp.h"
#include "capoBin.h"
#include "capoPlacer.h"
#include "partProbForCapo.h"
#include "splitLarge.h"
#include "splitSmall.h"
#include "subHGForCapo.h"

using std::cout;
using std::endl;
using std::max;
using std::min;
using std::vector;

PartitioningProblemForCapo::PartitioningProblemForCapo(
    const CapoBin& bin, const CapoPlacer& capo, BitBoard& edgesVisited,
    vector<vector<double> >& weights, BaseBinSplitter& splitter, bool isHCut,
    PartProbSizeType size, bool isRepart, const double splitPt,
    bool splitPtFixed, unsigned splitRow, double p0Capacity, Verbosity verb,
    bool honourGrpConstr)
    : PartitioningProblem(),
      _cellToBinMap(capo.getCellToBinMap()),
      _capoPl(capo),
      _size(size),
      _repartSmallWS(false),
      _splitPtFixed(splitPtFixed),
      _changedTol(false),
      _splitRow(splitRow),
      _xSplit(splitPt),
      _p0Capacity(p0Capacity),
      _edgesVisited(edgesVisited),
      _weights(weights),
      _verb(verb),
      _isRepart(isRepart) {
  abkfatal(
      size != Unknown,
      "This constructor must know the size classification of this problem");

  _bins = vector<const CapoBin*>(1);
  _bins[0] = &bin;

  initToNULL();

  // partitions get set only so we can construct the HGraph
  // and do termi-prop..it is expected that someone
  // will set the capacities, etc. before use
  const BBox& binBBox = bin.getBBox();
  _partitions = new std::vector<BBox>(2, BBox());
  if (isHCut) {
    _nDim = 2;
    _mDim = 1;
    (*_partitions)[0] = BBox(binBBox.xMin, splitPt, binBBox.xMax, binBBox.yMax);
    (*_partitions)[1] = BBox(binBBox.xMin, binBBox.yMin, binBBox.xMax, splitPt);
  } else {
    _nDim = 1;
    _mDim = 2;
    (*_partitions)[0] = BBox(binBBox.xMin, binBBox.yMin, splitPt, binBBox.yMax);
    (*_partitions)[1] = BBox(splitPt, binBBox.yMin, binBBox.xMax, binBBox.yMax);
  }

  _setPartitionCenters();

  // build the hgraph
  buildHGraphForPartitioning(honourGrpConstr);

  // do bloating if necessary
  double cellArea = bin.getTotalCellArea();

  if (_ownsHGraphs && capo.getParams().useBloating && size != Small) {
    if (capo.getParams().verb.getForActions() > 1)
      cout << " Cell Bloating " << cellArea << " -> ";
    cellBloating();
    cellArea = 0;
    for (vector<unsigned>::const_iterator nItr = bin.cellIdsBegin();
         nItr != bin.cellIdsEnd(); ++nItr) {
      unsigned originalCellID = (*nItr);

      const HGFNode& node = _subHG->getNewNodeByOrigIdx(originalCellID);
      cellArea += _subHG->getWeight(node.getIndex());
    }
    // icky, fix me
    CapoBin& thisBin = const_cast<CapoBin&>(bin);
    thisBin.setTotalCellArea(cellArea);
    if (capo.getParams().verb.getForActions() > 1) cout << cellArea << endl;
  }

  // set up capacities and tolerances
  _capacities =
      new std::vector<std::vector<double> >(2,
                                            std::vector<double>(1, 0.5 * cellArea));
  _totalWeight = new std::vector<double>(1, cellArea);
  _maxCapacities =
      new std::vector<std::vector<double> >(2, std::vector<double>(1, 0));
  _minCapacities =
      new std::vector<std::vector<double> >(2, std::vector<double>(1, 0));
  setTolerance(0);

  if (isHCut) {
    if (size == Large) {
      if (isRepart) {
        calculateCapacities_LargeH_repart();
        calculateTolerance_LargeH_repart();
      } else {
        calculateCapacities_LargeH();
        calculateTolerance_LargeH();
      }
    } else if (size == Medium) {
      if (isRepart) {
        calculateCapacities_SmallH();
        calculateTolerance_SmallH_repart(true);
      } else {
        calculateCapacities_SmallH();
        calculateTolerance_SmallH(true);
      }
    } else {
      if (isRepart) {
        calculateCapacities_SmallH();
        calculateTolerance_SmallH_repart(false);
      } else {
        calculateCapacities_SmallH();
        calculateTolerance_SmallH(false);
      }
    }
  } else {
    if (size == Large) {
      if (isRepart) {
        calculateCapacities_LargeV_repart();
        calculateTolerance_LargeV_repart();
      } else {
        calculateCapacities_LargeV();
        calculateTolerance_LargeV();
      }
    } else if (size == Medium) {
      if (isRepart) {
        calculateCapacities_SmallV();
        calculateTolerance_SmallV_repart(true);
      } else {
        calculateCapacities_SmallV();
        calculateTolerance_SmallV(true);
      }
    } else {
      if (isRepart) {
        calculateCapacities_SmallV();
        calculateTolerance_SmallV_repart(false);
      } else {
        calculateCapacities_SmallV();
        calculateTolerance_SmallV(false);
      }
    }
  }

  if (capo.getParams().splitterParams.useACG && size == Large && !isRepart &&
      bin.getWhitespace() >= 0. && capo.getMainLayerNum() <= 5) {
    const vector<unsigned>& cells_in_bin = bin.getCellIds();
    double accumX = 0.0, accumY = 0.0, accumArea = 0.0, minCellArea = DBL_MAX;
    // defensive driving:
    // if there are cells with zero area then all weights become 1
    for (unsigned i = 0; i < cells_in_bin.size(); ++i) {
      double currCellArea = capo.getNetlistHGraph().getNodeHeight(i) *
                            capo.getNetlistHGraph().getNodeWidth(i);
      if (minCellArea > currCellArea) {
        minCellArea = currCellArea;
      }
    }
    if (minCellArea == 0)
      abkwarn(
          0,
          "Applying ACG to bin with 0 area cells, overriding all areas to 1");

    for (unsigned i = 0; i < cells_in_bin.size(); ++i) {
      double currCellArea = capo.getNetlistHGraph().getNodeHeight(i) *
                            capo.getNetlistHGraph().getNodeWidth(i);
      // defensive driving:
      // if there are cells with zero area then all weights become 1
      if (minCellArea == 0) currCellArea = 1;
      Point currCellLoc = capo.getPlacement()[i];
      bin.getBBox().coerce(currCellLoc);
      if (!bin.getBBox().contains(currCellLoc)) {
        cout << "Cell not * in its bin BBox" << endl;
        cout << "\tBin BBox : " << bin.getBBox() << endl;
        cout << "\tCell Location : " << capo.getPlacement()[i] << endl;
      }

      accumX += currCellLoc.x * currCellArea;
      accumY += currCellLoc.y * currCellArea;
      accumArea += currCellArea;
    }
    double cogX = accumX / accumArea;
    double cogY = accumY / accumArea;
    cout << "Part prob for bin " << bin.getIndex() << " has COG (" << cogX
         << "," << cogY << ")" << endl;
    Point COG(cogX, cogY);
    if (!bin.getBBox().contains(COG)) {
      cout << "ACG PROBLEM!" << endl;
      cout << "Cell COG not in its bin BBox" << endl;
      cout << "\tBin BBox : " << bin.getBBox() << endl;
      cout << "\tCOG Location : " << COG << endl;
    }
    double AR = bin.getBBox().getWidth() / bin.getBBox().getHeight();
    double totCellArea = bin.getTotalCellArea();
    double H = sqrt(totCellArea / AR);
    double W = AR * H;
    BBox ACGrect(COG.x - 0.5 * W, COG.y - 0.5 * H, COG.x + 0.5 * W,
                 COG.y + 0.5 * H);
    cout << "ACGrect = " << ACGrect << endl;
    cout << "binBBox = " << bin.getBBox() << endl;

    // Coerce the ACG rect into the bin BBox
    double lDiff = bin.getBBox().xMin - ACGrect.xMin;
    if (lDiff > 0) {
      ACGrect.xMin += lDiff;
      ACGrect.xMax += lDiff;
    }
    double rDiff = ACGrect.xMax - bin.getBBox().xMax;
    if (rDiff > 0) {
      ACGrect.xMin -= rDiff;
      ACGrect.xMax -= rDiff;
    }
    double bDiff = bin.getBBox().yMin - ACGrect.yMin;
    if (bDiff > 0) {
      ACGrect.yMin += bDiff;
      ACGrect.yMax += bDiff;
    }
    double tDiff = ACGrect.yMax - bin.getBBox().yMax;
    if (tDiff > 0) {
      ACGrect.yMin -= tDiff;
      ACGrect.yMax -= tDiff;
    }
    cout << "ACGrect coerced = " << ACGrect << endl;
    cout << "splitPt = " << splitPt << endl;
    if (!bin.getBBox().contains(ACGrect)) {
      cout << "ACG PROBLEM! ACGrect is larger than bin bbox" << endl;
      cout << " binBBox: " << bin.getBBox() << endl;
      cout << " coerced: " << ACGrect << endl;
      abkwarn(0, "Problem in the computation of ACG 3, please report this bug");
      // Clean up the mess, figure out why some day
      if (ACGrect.getWidth() > bin.getBBox().getWidth()) {
        ACGrect.xMin = bin.getBBox().xMin;
        ACGrect.xMax = bin.getBBox().xMax;
      }
      if (ACGrect.getHeight() > binBBox.getHeight()) {
        ACGrect.yMin = bin.getBBox().yMin;
        ACGrect.yMax = bin.getBBox().yMax;
      }
    }

    double p0TargetPortion = 0.0;
    double p1TargetPortion = 0.0;
    if (isHCut) {
      p0TargetPortion = (ACGrect.yMax - splitPt) / ACGrect.getHeight();
      p1TargetPortion = (splitPt - ACGrect.yMin) / ACGrect.getHeight();
    } else {
      p0TargetPortion = (splitPt - ACGrect.xMin) / ACGrect.getWidth();
      p1TargetPortion = (ACGrect.xMax - splitPt) / ACGrect.getWidth();
    }
    cout << "Basic Target Portions 0: " << p0TargetPortion
         << " 1: " << p1TargetPortion << endl;
    p0TargetPortion = min(p0TargetPortion, 1.0);
    p0TargetPortion = max(p0TargetPortion, 0.0);
    p1TargetPortion = min(p1TargetPortion, 1.0);
    p1TargetPortion = max(p1TargetPortion, 0.0);
    cout << "Normalized Target Portions 0: " << p0TargetPortion
         << " 1: " << p1TargetPortion << endl;
    abkfatal(fabs(p0TargetPortion + p1TargetPortion - 1.0) < 1e-6,
             "ACG portions do not add to 1");

    cout << " Before ACG tolerance " << endl;
    cout << " Targets: " << getCapacities()[0][0] << "/"
         << getCapacities()[1][0] << endl;
    cout << " Max    : " << getMaxCapacities()[0][0] << "/"
         << getMaxCapacities()[1][0] << endl;
    cout << " Min    : " << getMinCapacities()[0][0] << "/"
         << getMinCapacities()[1][0] << endl;
    cout << " Tolerance : " << 100 * getTolerance(0) << endl;
    cout << endl;

    double p0Target = p0TargetPortion * totCellArea;
    double p1Target = p1TargetPortion * totCellArea;

    // if(p0Target > getMaxCapacities()[0][0] || p0Target <
    // getMinCapacities()[0][0])
    //{
    //   p0Target = min(p0Target, getMaxCapacities()[0][0]);
    //   p0Target = max(p0Target, getMinCapacities()[0][0]);
    //   p1Target = totCellArea - p0Target;
    // }

    // if(p1Target > getMaxCapacities()[1][0] || p1Target <
    // getMinCapacities()[1][0])
    //{
    //   p1Target = min(p1Target, getMaxCapacities()[1][0]);
    //   p1Target = max(p1Target, getMinCapacities()[1][0]);
    //   p0Target = totCellArea - p1Target;
    // }

    setCapacity(0, p0Target);
    setCapacity(1, p1Target);

    setTolerance(capo.getParams().splitterParams.defaultTolerance);

    cout << " After ACG tolerance " << endl;
    cout << " Targets: " << getCapacities()[0][0] << "/"
         << getCapacities()[1][0] << endl;
    cout << " Max    : " << getMaxCapacities()[0][0] << "/"
         << getMaxCapacities()[1][0] << endl;
    cout << " Min    : " << getMinCapacities()[0][0] << "/"
         << getMinCapacities()[1][0] << endl;
    cout << " Tolerance : " << 100 * getTolerance(0) << endl;
    cout << endl;
  }

  if (capo.getParams().splitterParams.WSOracle) {
    abkfatal(capo.getOraclePlacement().size() == capo.getPlacement().size(),
             "Oracle placement didn't get initialized");
    abkwarn(bin.getOracleCellIds().size() > 0,
            "Splitting a bin where the oracle has no cells");

    if (bin.getOracleCellIds().size() > 0) {
      double p0CellArea = 0., p1CellArea = 0.;
      if (isHCut) {
        for (vector<unsigned>::const_iterator i =
                 bin.getOracleCellIds().begin();
             i != bin.getOracleCellIds().end(); ++i) {
          if (capo.getOraclePlacement()[*i].y >= splitPt)
            p0CellArea += _capoPl.getNetlistHGraph().getWeight(*i);
          else
            p1CellArea += _capoPl.getNetlistHGraph().getWeight(*i);
        }
      } else {
        for (vector<unsigned>::const_iterator i =
                 bin.getOracleCellIds().begin();
             i != bin.getOracleCellIds().end(); ++i) {
          if (capo.getOraclePlacement()[*i].x >= splitPt)
            p1CellArea += _capoPl.getNetlistHGraph().getWeight(*i);
          else
            p0CellArea += _capoPl.getNetlistHGraph().getWeight(*i);
        }
      }

      double p0Target =
          (bin.getTotalCellArea() * p0CellArea) / (p0CellArea + p1CellArea);
      double p1Target = bin.getTotalCellArea() - p0Target;

      setCapacity(0, p0Target);
      setCapacity(1, p1Target);

      setTolerance(capo.getParams().splitterParams.WSOracleTol);
    }
  }

  double curTol = getTolerance();
  if (curTol > _capoPl.getParams().splitterParams.tolCap) {
    multiplyTolerance(_capoPl.getParams().splitterParams.tolCap / curTol);
  }

  // added by jflu & aaronnn
  if (_capoPl.getParams().prePartFix && size != Small && !isRepart) {
    fixCellsUsingIntermediatePlacement();
  }

  setBuffer();

  _ownsData = false;
}

PartitioningProblemForCapo::PartitioningProblemForCapo(
    const vector<const CapoBin*>& bins, const CapoPlacer& capo,
    const vector<const CapoBin*>& cellToBinMap, BitBoard& edgesVisited,
    vector<vector<double> >& weights, Verbosity verb)
    : PartitioningProblem(),
      _bins(bins),
      _cellToBinMap(cellToBinMap),
      _capoPl(capo),
      _size(Unknown),
      _repartSmallWS(false),
      _splitPtFixed(false),
      _changedTol(false),
      _splitRow(UINT_MAX),
      _xSplit(0.),
      _p0Capacity(0.),
      _edgesVisited(edgesVisited),
      _weights(weights),
      _hCutWSMethod(UINT_MAX),
      _optimisticAdj(UINT_MAX),
      _verb(verb),
      _isRepart(false) {
  abkfatal(_bins.size() == 2, "can only build problems for 2 bins");

  initToNULL();

  setPartDims();  // use the bin shapes/capacities for the partitions

  buildUnweightedHGraph(false);

  setBuffer();

  _ownsData = false;
}

PartitioningProblemForCapo::PartitioningProblemForCapo(
    const vector<const CapoBin*>& bins, const CapoPlacer& capo,
    const vector<const CapoBin*>& cellToBinMap, BitBoard& edgesVisited,
    vector<vector<double> >& weights, const vector<double>& partCaps,
    Verbosity verb)
    : PartitioningProblem(),
      _bins(bins),
      _cellToBinMap(cellToBinMap),
      _capoPl(capo),
      _size(Unknown),
      _repartSmallWS(false),
      _splitPtFixed(false),
      _changedTol(false),
      _splitRow(UINT_MAX),
      _xSplit(0.),
      _p0Capacity(0.),
      _edgesVisited(edgesVisited),
      _weights(weights),
      _hCutWSMethod(UINT_MAX),
      _optimisticAdj(UINT_MAX),
      _verb(verb),
      _isRepart(false) {
  abkfatal(_bins.size() == 2, "can only build problems for 2 bins");

  initToNULL();

  setPartDims(partCaps);  // use the bin shapes/capacities for the partitions

  buildUnweightedHGraph(false);

  setBuffer();

  _ownsData = false;
}

PartitioningProblemForCapo::PartitioningProblemForCapo(
    PartitioningProblemForCapo& origPartProb, bool shred)
    : PartitioningProblem(),
      _cellToBinMap(origPartProb._capoPl.getCellToBinMap()),
      _capoPl(origPartProb._capoPl),
      _size(origPartProb._size),
      _repartSmallWS(origPartProb._repartSmallWS),
      _splitPtFixed(origPartProb._splitPtFixed),
      _changedTol(origPartProb._changedTol),
      _splitRow(origPartProb._splitRow),
      _xSplit(origPartProb._xSplit),
      _p0Capacity(origPartProb._p0Capacity),
      _edgesVisited(origPartProb._edgesVisited),
      _weights(origPartProb._weights),
      _verb(origPartProb._verb),
      _isRepart(origPartProb._isRepart) {
  _bins = vector<const CapoBin*>(1);
  _bins[0] = origPartProb._bins[0];
  initToNULL();
  double cellArea = _bins[0]->getTotalCellArea();

  _capacities =
      new std::vector<std::vector<double> >(
          2, std::vector<double>(1, cellArea / 2.0));
  _totalWeight = new std::vector<double>(1, cellArea);
  _maxCapacities =
      new std::vector<std::vector<double> >(2, std::vector<double>(1, 0));
  _minCapacities =
      new std::vector<std::vector<double> >(2, std::vector<double>(1, 0));
  const std::vector<std::vector<double> >& origCap =
      origPartProb.getCapacities();
  const std::vector<std::vector<double> >& origMaxCap =
      origPartProb.getMaxCapacities();
  const std::vector<std::vector<double> >& origMinCap =
      origPartProb.getMinCapacities();
  for (unsigned i = 0; i < (*_capacities).size(); ++i)
    for (unsigned j = 0; j < (*_capacities)[i].size(); ++j) {
      (*_capacities)[i][j] = origCap[i][j];
      (*_maxCapacities)[i][j] = origMaxCap[i][j];
      (*_minCapacities)[i][j] = origMinCap[i][j];
    }
  _partitions = new std::vector<BBox>(2, BBox());
  const std::vector<BBox>& origPartitions = origPartProb.getPartitions();
  for (unsigned i = 0; i < origPartitions.size(); ++i)
    setPartition(i, origPartitions[i]);
  _nDim = origPartProb._nDim;
  _mDim = origPartProb._mDim;
  _setPartitionCenters();

  buildShredHGraph(origPartProb);

  _ownsData = origPartProb.isDataOwned();
}

PartitioningProblemForCapo::~PartitioningProblemForCapo() {
  destroyHGraph();

  if (_fixedConstr != NULL) {
    delete _fixedConstr;
    _fixedConstr = NULL;
  }

  if (_solnBuffers != NULL) {
    delete _solnBuffers;
    _solnBuffers = NULL;
  }

  if (_capacities != NULL) {
    delete _capacities;
    _capacities = NULL;
  }

  if (_partitions != NULL) {
    delete _partitions;
    _partitions = NULL;
  }

  if (_padBlocks != NULL) {
    delete _padBlocks;
    _padBlocks = NULL;
  }
  if (_terminalToBlock != NULL) {
    delete _terminalToBlock;
    _terminalToBlock = NULL;
  }
}
