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

#ifndef _BASEBINSPLITTER_H_
#define _BASEBINSPLITTER_H_

#include <cfloat>
#include <climits>
#include <utility>

#include <string>
#include <vector>
#include "CongestionMaps/ISPD04CongMap.h"
#include "capoBin.h"
#include "capoParams.h"
#include "capoPlacer.h"
#include "partProbForCapo.h"

class BaseBinSplitter {
 protected:
  CapoBin& _bin;
  std::vector<CapoBin*> _newBins;
  unsigned _splitRow;
  double _xSplit;
  CapoPlacer& _capo;

  CapoSplitterParams _params;

 public:
  unsigned findBestSplitRow(const std::vector<double>& cellAreas,
                            double p0WSWeight, double p1WSWeight,
                            std::vector<double>& actualCaps) const;

 protected:
  unsigned shiftSplitRow(unsigned curRow, const std::vector<double>& cellAreas,
                         std::vector<double>& actualCaps,
                         double localWS) const;

  // <aaronnn> adjust WS of split depending in utilization density
  // (e.g. one large macro is 'denser' than small cells with similar total area)
  double adjustSplitByUtilizationDensity(
      const PartitioningProblemForCapo& problem, double splitPt, bool isHCut);

  // shiftSplitRowAndSnap(...) is a shell function for calling
  // shiftSplitRow(...) followed by a call to snapRowToOracleFeature(...)
  // if useCutOracle is not enabled, then shiftSplitRowAndSnap will accept
  // the result of shiftSplitRow(...) and not snap
  unsigned shiftSplitRowAndSnap(unsigned curRow,
                                const std::vector<double>& cellAreas,
                                std::vector<double>& actualCaps,
                                double localWS, double roundPct) const;

  unsigned snapRowToOracleFeature(unsigned curRow,
                                  const std::vector<double>& cellAreas,
                                  std::vector<double>& actualCaps,
                                  double localWS, double roundPct) const;

  // create bins after a horizontal split
  // first version finds the best splitRow..the
  // second uses the one it's given
  void createHSplitBins(bool isLarge, PartitioningProblemForCapo& problem);
  void createHSplitBins(bool isLarge, PartitioningProblemForCapo& problem,
                        unsigned splitRow);

 public:
  void createHSplitBins(
      const std::vector<std::vector<unsigned> >& cellsInNewBins,
      unsigned splitRow, double p0MaxCap, double p1MaxCap);

 protected:
  // create bins after a vertical split
  void createVSplitBins(bool isLarge, PartitioningProblemForCapo& problem,
                        double xSplit);

 public:
  void createVSplitBins(
      const std::vector<std::vector<unsigned> >& cellsInNewBins,
      double xSplit, double p0MaxCap, double p1MaxCap);

  BaseBinSplitter(CapoBin& binToSplit, CapoPlacer& capo,
                  CapoSplitterParams params)
      : _bin(binToSplit),
        _splitRow(UINT_MAX),
        _xSplit(DBL_MAX),
        _capo(capo),
        _params(params) {}

  virtual ~BaseBinSplitter() {}  // does not own new bins

  std::vector<CapoBin*>& getNewBins() { return _newBins; }

  unsigned getSplitRow() const { return _splitRow; }
  double getXSplit() const { return _xSplit; }

  // give the partition areas in presence of filler cells
  bool partAreasForFillerCells(bool splitDir, std::vector<double>& partAreas);

  std::pair<double, double> getCongestionRatiosH(
      unsigned splitRow, const ISPD04CongMap& congmap) const;
  std::pair<double, double> getCongestionRatiosV(
      double xSplit, const ISPD04CongMap& congmap) const;

  bool shouldCongestionShiftH(unsigned splitRow, const ISPD04CongMap& congmap,
                              double avgCong) const;
  bool shouldCongestionShiftV(double xSplit, const ISPD04CongMap& congmap,
                              double avgCong) const;

  void shiftOraclePlacementH(double ecoSplitRowCoord,
                             double shiftedSplitRowCoord);
  void shiftOraclePlacementV(double ecoXSplit, double shiftedXSplit);

  double partitionByGroups(PartitioningProblemForCapo& problem);
  bool shouldPartitionByGroups(const PartitioningProblemForCapo& problem) const;
};

#endif
