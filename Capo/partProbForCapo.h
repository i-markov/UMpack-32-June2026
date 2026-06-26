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

#ifndef __CAPO_PART_PROB_H__
#define __CAPO_PART_PROB_H__

#include <iostream>

#include "Partitioning/partProb.h"
#include "baseBinSplitter.h"
#include "subHGForCapo.h"

class BaseBinSplitter;
class Partitioning;
class CapoBin;
class SubHGraph;

enum PartProbSizeType { Small, Medium, Large, Unknown };

class PartitioningProblemForCapo : public PartitioningProblem {
  friend class CapoPlacer;

 private:
  // utility function added by DAP to simplify implementation of the
  // of weighted terminal propagation algorithm
  void addEdgePropagateTerminalsReweight(
      SubHGraphForCapo* newHG, const itHGFEdgeLocal& eItr,
      const std::vector<unsigned>& nonTermsOnNet, bool term0, bool term1,
      double reweightingFactor);

 protected:
  std::vector<const CapoBin*> _bins;
  const std::vector<const CapoBin*>& _cellToBinMap;
  const CapoPlacer& _capoPl;
  PartProbSizeType _size;
  unsigned _nDim;
  unsigned _mDim;
  bool _shouldBeRepartitioned;
  bool _repartSmallWS;
  bool _splitPtFixed;
  bool _changedTol;
  unsigned _splitRow;
  double _xSplit;
  double _p0Capacity;

  std::vector<int> _yRows;

  SubHGraph* _subHG;  // IMPORTANT: this points to the
                      // same object as the _hgraph value
                      // inherited from the base. But, it's
                      //'subHG' properties are necessary, so this
                      // allows for using them.

  BitBoard& _edgesVisited;
  std::vector<std::vector<double> >& _weights;

  void setTolerances();
  void setBuffer();
  void setPartDims();
  void setPartDims(const std::vector<double>& partCaps);

  void buildHGraphForPartitioning(bool honourGrpConstr = false);
  void computeUnifiedTheto();
  void constructHG(bool honourGrpConstr);
  void cellBloating();

  void buildUnweightedHGraph(bool honourGrpConstr);
  // following ctor shreds the big macros. assumes tolerances are set aprior
  void buildShredHGraph(PartitioningProblemForCapo& origPartProb);

  unsigned _hCutWSMethod;
  double _optimisticAdj;

  Verbosity _verb;

  bool _ownsHGraphs;
  bool _isRepart;

 public:
  // simpler method used by the 'SplitSmallBin' methods.
  // does not set tolerances, etc.
  PartitioningProblemForCapo(const CapoBin& bin, const CapoPlacer& capo,
                             BitBoard& edgesVisited,
                             std::vector<std::vector<double> >& weights,
                             BaseBinSplitter& splitter, bool isHCut,
                             PartProbSizeType size, bool isRepart,
                             const double splitPt, bool splitPtFixed,
                             unsigned splitRow, double p0Capacity,
                             Verbosity verb, bool honourGrpConstr);

  // build a problem to repartition 2 bins
  PartitioningProblemForCapo(const std::vector<const CapoBin*>& bins,
                             const CapoPlacer& capo,
                             const std::vector<const CapoBin*>& cellToBinMap,
                             BitBoard& edgesVisited,
                             std::vector<std::vector<double> >& weights,
                             Verbosity verb = Verbosity("1 1 1"));

  // build a problem to repartition 2 bins given congestion info
  PartitioningProblemForCapo(const std::vector<const CapoBin*>& bins,
                             const CapoPlacer& capo,
                             const std::vector<const CapoBin*>& cellToBinMap,
                             BitBoard& edgesVisited,
                             std::vector<std::vector<double> >& weights,
                             const std::vector<double>& partCaps,
                             Verbosity verb = Verbosity("1 1 1"));

  // this ctor builds copy of a partitioning problem and
  // shreds the big macros into a grid
  PartitioningProblemForCapo(PartitioningProblemForCapo& origPartProb,
                             bool shred = true);

  virtual ~PartitioningProblemForCapo();

  unsigned getNDim() const { return _nDim; }
  unsigned getMDim() const { return _mDim; }

  SubHGraph& getSubHGraph() const { return *_subHG; }

  bool getRepartSmallWS() const { return _repartSmallWS; }
  bool getSplitPtFixed() const { return _splitPtFixed; }
  bool getChangedTol() const { return _changedTol; }
  unsigned getSplitRow() const { return _splitRow; }
  void setSplitRow(unsigned r) { _splitRow = r; }
  double getXSplit() const { return _xSplit; }
  void setXSplit(double x) { _xSplit = x; }
  double getP0Capacity() const { return _p0Capacity; }
  void setP0Capacity(double c) { _p0Capacity = c; }

  void setPartition(unsigned partNum, const BBox& newBBox) {
    (*_partitions)[partNum] = newBBox;
  }

  void setCapacity(unsigned partNum, double newCap) {
    (*_capacities)[partNum][0] = newCap;
  }

  void setMaxCapacity(unsigned partNum, double newCap) {
    (*_maxCapacities)[partNum][0] = newCap;
  }

  void setMinCapacity(unsigned partNum, double newCap) {
    (*_minCapacities)[partNum][0] = newCap;
  }

  // resetTolerance is to be used when partitioning a single
  // bin, after the initial partitioning
  void resetTolerance(double s0, double s1);

  bool shouldBeRepartitioned() const { return _shouldBeRepartitioned; }

  void setTolerance(double tol);
  void increaseTolerance(double tol);
  void decreaseTolerance(double tol);
  void multiplyTolerance(double factor);
  double getTolerance(unsigned wtNum = 0) const;
  void setSafeWSToleranceH(unsigned splitRow, double safeWS);
  void setMinLocalWSToleranceH(unsigned splitRow, double minLocalWS);
  void setPaperMethodToleranceH(unsigned splitRow);
  void setSafeWSToleranceV(double xSplit, double safeWS);
  void setMinLocalWSToleranceV(double xSplit, double minLocalWS);
  void setPaperMethodToleranceV(double xSplit);

  // jflu

  // Small Horizontal
  void calculateCapacities_SmallH();
  void calculateTolerance_SmallH(bool midSize);
  void calculateTolerance_SmallH_repart(bool midSize);
  // Small Vertical
  void calculateCapacities_SmallV();
  void calculateTolerance_SmallV(bool midSize);
  void calculateTolerance_SmallV_repart(bool midSize);
  // Large Horizontal
  void calculateCapacities_LargeH();
  void calculateCapacities_LargeH_repart();
  void calculateTolerance_LargeH();
  void calculateTolerance_LargeH_repart();
  // Large Vertical
  void calculateCapacities_LargeV();
  void calculateCapacities_LargeV_repart();
  void calculateTolerance_LargeV();
  void calculateTolerance_LargeV_repart();

  bool hPartCapsFrmRgnConstr(unsigned splitRow, std::vector<double>& partCaps);
  bool vPartCapsFrmRgnConstr(double xSplit, std::vector<double>& partCaps);

  // jflu and aaronnn
  void fixCellsUsingIntermediatePlacement();

  void destroyHGraph();
  void swapHGraphs(PartitioningProblemForCapo& newProb);
};

#endif
