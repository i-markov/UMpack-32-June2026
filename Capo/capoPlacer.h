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

// Created: Sep 11, 1997 by Igor Markov & Andy Caldwell

// CHANGES
// 971113 aec added class CapoPlacer::Parameters
// 971118 aec added SA multi-start and multi-level to Parameters
// 980313 ilm added CapoPlacer::Parameters ctor from argc, argv
//            reworked CapoPlacer::Parameters and its output
//            now the original argc, argv are saved in the class
// 980325 ilm added n,m and partFuzziness to CapoPlacer::Parameters
// 980401 ilm added savePartProb parameter
// 980401 ilm added capoPlacer::splitBin()
// 980402 aec added capoPlacer::runSAPlace()
// 980609 ilm changed ClusterTree to ClusteringFromDB
// 980619 aec updated to ClusterTreeFromDB
// 990317 ilm bin is not passed as const to constructNewBins
//	     cause its neighbors are unlinked

#ifndef __CAPOPLACER_H__
#define __CAPOPLACER_H__

#include <cfloat>
#include <climits>
#include <cstddef>
#include <functional>
#include <iosfwd>
#include <unordered_map>
#include <utility>

#include "ABKCommon/abkassert.h"
#include "ABKCommon/abkrand.h"
#include "ABKCommon/abktypes.h"
#include <string>
#include <vector>
#include "Constraints/constraints.h"
#include "Ctainers/bitBoard.h"
#include "HGraph/subHGraph.h"
#include "HGraphWDims/hgWDims.h"
#include "MLPart/mlPart.h"
#include "RBPlace/RBPlacement.h"
#include "SmallPlacers/baseSmallPl.h"
#include "baseBinSplitter.h"
#include "capoBin.h"
#include "capoParams.h"
// #include "RBPlace/slicingTree.h"
#include "CongestionMaps/CongestionMaps.h"
#include "CongestionMaps/ISPD04CongMap.h"
#include "GeomTrees/FastSteiner/global.h"
#include "Geoms/ptsetpt.h"

class PartitioningProblemForCapo;
class CapoSplitterParams;

struct _splitData {
  bool splitHoriz;
  unsigned splitRow;
  double xSplit;
  CapoBin *parent, *p0, *p1;
};
typedef struct _splitData splitData;

class CapoPlacer {
  friend class LookAheadSplitter;
  friend class SplitRowOfBins;

 private:
  CapoParameters _params;

  RBPlacement& _rbplace;

  BBox _coreBBox;
  PlacementWOrient _placement;
  // Capo does not populate the
  // RBPlacement incrementally..but
  // rather, copies this placement into
  // the RBPlace all at once at the end
  PlacementWOrient _oraclePlacement;

  // added by sadya to include CongestionMaps
  CongestionMaps* _congestionMap;
  // added by royj for a different congestion map
  ISPD04CongMap* _ispd04CongMap;
  PlacementWOrient* _ispd04CongMapPlacement;

  bit_vector _placed;  // has the node been assigned a location

  std::vector<CapoBin*> _layout[2];
  std::vector<CapoBin*>* _curLayout;
  std::vector<CapoBin*> _placedBins;
  std::vector<CapoBin*> _saveBins;

  std::vector<splitData> _splits;

  unsigned _totalLayerNum, _mainLayerNum, _feedbackLayerNum;

 public:
  unsigned nextBinNum;

 private:
  std::vector<const CapoBin*> _cellToBinMap;
  // contains a pointer to the bin containing each cell.
  // cells not in a bin have the pointer == NULL

  std::vector<BBox> _netUpperBounds;
  std::vector<BBox> _netLowerBounds;

  // useful bit-vectors (used when traversing nets/edges).
  // users must 'clear' the vectors before use
  bit_vector _nodeSeen;
  bit_vector _edgeSeen;

  // <aaronnn> mark nodes needing FP to fix after FP
  bit_vector _nodesNeedingFP;

 public:
  // royj: bitboards used by small placement and partitioning
  // that were once static but are now just owned by Capo
  SmallPlacementBitBoardContainer spBitBoards;
  BBPartBitBoardContainer bbBitBoards;

 private:
  // jflu: pointSet to improve net traversal speed
  std::vector<std::vector<PtSetPoint> > _pointSet_fixed;
  std::vector<std::vector<PtSetPoint> > _pointSet_movable;

 public:
  BitBoard partProbEdgesVisited;
  std::vector<std::vector<double> > thetoWeights;

 private:
  // Data tracking structures

  std::vector<double> _estBBoxWLPerLayer;  // est BBox WL at each level
  std::vector<unsigned> _totCutPerLayer;   // total cut at each level

  std::vector<unsigned> _nodesInEachBin;  // used for calculating rent #s
  std::vector<unsigned> _terminalsToEachBin;

  std::vector<unsigned> _numOrigNets;
  std::vector<unsigned> _numEssentialNets;
  std::vector<unsigned> _origNetPins;
  std::vector<unsigned> _essentialNetPins;
  std::vector<unsigned> _numProblemsOfSize;

  //'cut density' data storage
  //'Stats' are #cut/length of cut
  std::vector<std::vector<std::pair<unsigned, double> > > _vCutStats;
  std::vector<std::vector<std::pair<unsigned, double> > > _hCutStats;
  //'Ratios' are #newly cut/num cut
  std::vector<std::vector<std::pair<unsigned, unsigned> > > _vCutRatios;
  std::vector<std::vector<std::pair<unsigned, unsigned> > > _hCutRatios;

  // Storage for statistics for cut nets
  std::unordered_map<unsigned, std::vector<unsigned>, std::hash<unsigned>,
           std::equal_to<unsigned> >
      _layerCut;
  bit_vector _netHasBeenCut;

  unsigned _numNotSolved;  // number of smallPart problems
  // not solved exactly
  unsigned _numSmallPartProbs;

  double _totalLayoutArea;
  double _totalOverfill;
  double _maxOverfill;

  class FloorplanStats {
   public:
    // updated in capoMainLoop.cxx
    static unsigned numInstances;
    static unsigned numSuccessfulInstances;
    static unsigned numFailedInstances;
    static unsigned instancesMaxMacroNum;
    static unsigned failedInstancesMaxMacroNum;
    static unsigned totalFPcells;
    static unsigned totalFPmacros;
    static unsigned totalFPobstacles;
    static std::vector<unsigned> cellsInInstance;
    static std::vector<unsigned> macrosInInstance;
    static std::vector<unsigned> obstaclesInInstance;
    static void clearFPStats() {
      numInstances = 0;
      numSuccessfulInstances = 0;
      numFailedInstances = 0;
      instancesMaxMacroNum = 0;
      failedInstancesMaxMacroNum = 0;
      totalFPcells = 0;
      totalFPmacros = 0;
      totalFPobstacles = 0;
      cellsInInstance.clear();
      macrosInInstance.clear();
      obstaclesInInstance.clear();
    }
  };

  RandomRawUnsigned SeedGen;

  enum LayerType { FP_LAYER, PATOMA_LAYER };

  // following stores the input desired slicing tree for the placement
  // In SYNP
  std::vector<std::pair<std::string, double> > _requiredSlicingTree;

  //  RBPlace::SlicingTree _slicingTree;

  //  const std::vector<char*>* _altCellNames;  // Capo does not own this
  //  pointer

  std::vector<BBox> fixedObstacles;
  std::vector<unsigned> relevantObstacles;
  std::vector<std::vector<std::vector<unsigned> > > obstacleGrid;

 public:
  class capoTimer {
   public:
    static double FMTime;
    static double MLTime;
#ifdef USEHMETIS
    static double HMetisTime;
#endif
    static double BBTime;
    static double ECOTime;
    static double SmTime;
    static double SmPlaceProb;
    static double SetupTime;
    static double FloorplanTime;
    static double FloorClusterTime;
    static double FloorAnnealTime;
    static double MacroOverlapTime;
    static double WLCalcTime;
    static double FeedbackProcessing;
    static double ThunderTime;
    static double SORTime;
    static double PointSetTime;
    static double FastPlaceMoveTime;
    static double BloatTime;
    static double CongMapTime;
    static double HGraphRebuildTime;
#ifdef USEFLOW
    static double FlowTime;
#endif
    static void clearTimer(void) {
      FMTime = 0.;
      MLTime = 0.;
#ifdef USEHMETIS
      HMetisTime = 0.;
#endif
      BBTime = 0.;
      ECOTime = 0.;
      SmTime = 0.;
      SmPlaceProb = 0.;
      SetupTime = 0.;
      FloorplanTime = 0.;
      FloorClusterTime = 0.;
      FloorAnnealTime = 0.;
      MacroOverlapTime = 0.;
      WLCalcTime = 0.;
      FeedbackProcessing = 0.;
      ThunderTime = 0.;
      SORTime = 0.;
      PointSetTime = 0.;
      FastPlaceMoveTime = 0.;
      BloatTime = 0.;
      CongMapTime = 0.;
      HGraphRebuildTime = 0.;
#ifdef USEFLOW
      FlowTime = 0.;
#endif
    }
  } capoTimes;

  class EcoStats {
   public:
    static unsigned ecoPartsUsed;
    static unsigned ecoPartsRejected;
    static unsigned ecoFPsUsed;
    static void clearEcoStats() {
      ecoPartsUsed = 0;
      ecoPartsRejected = 0;
      ecoFPsUsed = 0;
    }
  };

  CapoPlacer(RBPlacement& rbplace, const CapoParameters& params,
             CongestionMaps* congestionMap = NULL);

  CapoPlacer(RBPlacement& rbplace, const char* hclFileName,
             const char* binPLFileName, const CapoParameters& params);

  ~CapoPlacer();

  double estimateWL();

  const std::vector<double>& getBBoxPerLayer() const {
    return _estBBoxWLPerLayer;
  }
  const std::vector<unsigned>& getCutPerLayer() const {
    return _totCutPerLayer;
  }

  const std::vector<unsigned>& getNodesInEachBin() const {
    return _nodesInEachBin;
  }
  const std::vector<unsigned>& getTerminalsToEachBin() const {
    return _terminalsToEachBin;
  }
  const std::vector<const CapoBin*>& getCellToBinMap() const {
    return _cellToBinMap;
  }

  BBox getPinOffset(unsigned cellId, unsigned netId) const;
  BBox getPinLocation(unsigned cellId, unsigned netId) const;

  const HGraphWDimensions& getNetlistHGraph() const {
    return _rbplace.getHGraph();
  }

  const CapoParameters& getParams() const { return _params; }

  const PlacementWOrient& getPlacement() const { return _placement; }
  const PlacementWOrient& getOraclePlacement() const {
    return _oraclePlacement;
  }
  PlacementWOrient& getMutableOraclePlacement() { return _oraclePlacement; }

  //  const std::vector<char*>* getAltCellNames() const { return _altCellNames;
  //  }

  const RBPlacement& getRBP() const { return _rbplace; }

  void destroyHGraph() const { _rbplace.destroyHGraph(); }

  void rebuildHGraph(void) const {
    Timer rebuildTime;
    bool rebuilt = _rbplace.rebuildHGraph();
    rebuildTime.stop();
    CapoPlacer::capoTimer::HGraphRebuildTime += rebuildTime.getUserTime();
    if (rebuilt && _params.verb.getForActions() > 0)
      std::cout << "Rebuilding the HGraph took " << rebuildTime << std::endl;
  }

  void clearNamesFromHGraph(void) const { _rbplace.clearNamesFromHGraph(); }

  //  RBPlace::SlicingTree getSlicingTree(void){ return _slicingTree; }

  unsigned getMainLayerNum(void) const { return _mainLayerNum; }

  // by sadya in SYNP
  const Constraints& constraints;  // referenced from constraints in RBlacement
  std::vector<std::pair<BBox, double> >&
      regionUtilization;  // by sadya for region
  // utilization constraints. added while at SYNP
  std::vector<std::pair<BBox, std::vector<unsigned> > >& groupRegionConstr;
  std::vector<unsigned>& cellToGrpMapping;  // do we really need this?? TO DO
  bool isCutlinePredetermined(const std::string name, double& loc,
                              bool& direction);
  void readCutlineFile(const std::string cutlineFileName);

  // <aaronnn> track nodes needing FP, to fix after FP
  bool nodeNeedsFP(unsigned nodeID) const;
  void markNodeNeedsFP(unsigned nodeID);
  void unmarkNodeNeedsFP(unsigned nodeID);
  void markNodesNeedFP(const std::vector<unsigned>& cellIds);
  void markBinNeedsFP(unsigned binID);

  // jflu: pointSet to speedup net traversal
  void buildPointSetsFromScratch();
  void updatePointSetsForMovedCell(unsigned cellId, const Point& oldPos,
                                   const Point& newPos);
  void updatePointSetsForNewlyPlacedCell(unsigned cellId, const Point& oldPos,
                                         const Point& newPos);

  const std::vector<PtSetPoint>& getPointSet_fixed(unsigned netID) const;
  const std::vector<PtSetPoint>& getPointSet_movable(unsigned netID) const;

  bool isPlaced(unsigned nodeID) const { return _placed[nodeID]; }

  // jflu: congestion estimation vector based on connectivity
  std::vector<double> _congestion;
  double _congestionCutoff;

  // jflu: congestion estimation based on connectivity
  void updateCongestionEstimation();

  const ISPD04CongMap* getISPD04CongMap(void) const { return _ispd04CongMap; }

 private:
  void setupAndCheck();
  // sets up the _coreCells
  // checks that Pads are not inside, and core cells are
  // not outside, the coreBBox.
  // Creates the first(top-level) CapoBin, and places it
  // in _layout[0]. sets _curLayout to _layout[0].

  void readBinsFile(const char* binFileName, const char* hclFileName);

  // <aaronnn> figure out if we're at the bottom layers
  bool atBottomLayers();

  // <aaronnn> perform lookAheadFP and report success/failure
  bool lookAheadFPSuccess(unsigned b, std::vector<unsigned> macrosInBin,
                          std::vector<CapoBin*>& curLayout,
                          std::vector<CapoBin*>& newLayout);

  enum AllowedPartDir { HOnly, VOnly, HAndV };

  // one layer of bin refinment. Returns true if any
  // bin remains unplaced
  bool doOneLayer(AllowedPartDir partDir, bool& fpdSomething);
  void doVOnlyLayer(std::vector<CapoBin*>& curLayout,
                    std::vector<CapoBin*>& newLayout);
  void doHOnlyLayer(std::vector<CapoBin*>& curLayout,
                    std::vector<CapoBin*>& newLayout);
  bool doHAndVLayer(std::vector<CapoBin*>& curLayout,
                    std::vector<CapoBin*>& newLayout);

  // does one layer of only FP runs during mixed-size placement
  // returns true if atleast one bin is FPed
  bool doFPLayer(std::vector<CapoBin*>& curLayout,
                 std::vector<CapoBin*>& newLayout);

  // does one layer of only FP runs using PATOMA instead of the
  // standard approach
#ifdef USEPATOMA
  bool doPatomaLayer(std::vector<CapoBin*>& curLayout,
                     std::vector<CapoBin*>& newLayout);
#endif
  // <aaronnn> do some calculations to decide which way to split
  bool getDefaultSplitDir(CapoBin& bin);

  // refineBin returns true iff bin was partitioned,
  // and its children appended to newLayout
  bool refineBin(CapoBin& bin, std::vector<CapoBin*>& newLayout,
                 /* unsigned layerNum, */ AllowedPartDir partDir,
                 unsigned minNumCellsToPart, bool useCongestionInfo = false);

  void doOverlapping(std::vector<CapoBin*>& layout);

  void plotBins(/* unsigned layerNum, */ std::vector<CapoBin*>& layout);
  // saves a gnuplot script capoBins-layerX.gp

  void saveFloorplanProblem(const char* baseFileName);
  // saves all capoBins in the layout as a floorplanning instance.
  // also saves a clustering file, which acts as a mapping between
  // the netlist cells, and the bins.

  Point getWeberLocation(const CapoBin& bin);
  Point getWeberLocation(unsigned nodeIdx);
  BBox getWeberRegion(const CapoBin& bin);
  BBox getWeberRegion(unsigned nodeIdx);
  unsigned nodesToMove(const CapoBin& bin);

  void preSeedPartitioning(const CapoBin& bin,
                           PartitioningProblemForCapo& probToSeed);
  // NOTE: preSeedPart has a side-effect.  It changes the partitionings
  // in probToSeed

  // sets _cellToBinMap and placement
  void updateInfoAboutBin(const CapoBin* bin);

  // only update placement info
  void updatePlInfoAboutBin(const CapoBin* bin);

  // update only mapping info, not placement. used for capoNE
  void updateMapInfoAboutBin(const CapoBin* bin);

  // calculate the nets cut this layer, and those cut before
  void collectNetStats(std::vector<HGFEdge*>& netsCutThisLayer,
                       std::vector<HGFEdge*>& netsCutBefore);

  void replaceSmallBins();

  void repartitionBins(CapoBin* bin0, CapoBin* bin1);

  void repartitionBins(CapoBin* bin0, CapoBin* bin1, double relCongestionBin0,
                       double relCongestionBin1);

  // floorplan this bin using the standard flow
  bool floorplanBin(CapoBin& bin, bool minWL = true, bool lookAheadFP = false);

  // floorplan bin with SOR
  bool floorplanBinWithSOR(CapoBin& bin);

#ifdef USEFLOW
  bool floorplanBinWithFlows(CapoBin& bin);
#endif

  CapoBin* mergeBinWithSibling(CapoBin& bin);

  // floorplan this bin using PATOMA
#ifdef USEPATOMA
  bool patomaBin(CapoBin& bin, bool minWL = true);
#endif

  // process the placed macros after every layer
  void processMacrosAfterlayer(std::vector<CapoBin*>& layout,
                               LayerType layerType = FP_LAYER);
  void getFPedBinsIntoNewLayout(std::vector<CapoBin*>& layout,
                                std::vector<CapoBin*>& newLayout);
#ifdef USEPATOMA
  void getPatomaedBinsIntoNewLayout(std::vector<CapoBin*>& layout,
                                    std::vector<CapoBin*>& newLayout);
#endif

  // get the bins with unplaced macros in front of the list
  void getBinsWMacrosInFront(std::vector<CapoBin*>& layout);

  // by royj for new procedures after Capo
  // void doOptimalInterleaving(void);
  bool doRedistribution(void);
  void swapCells(CapoBin* b1, unsigned b1Cell, CapoBin* b2, unsigned b2Cell);
  CapoBin* findNearestBin(std::vector<CapoBin*> haveSpace,
                          const CapoBin& start, double spaceneeded,
                          unsigned* swapCell, double* bestdistance);

 public:
  void printCutDensityStats();
  void printSmallProblemStats();
  void printNetlistStats();

  void getNetCutInfo(std::vector<unsigned>& externalCutNets,
                     std::vector<unsigned>& containedNets) const;
  unsigned getTotalNetCut() const;

  // Begin functions added by DP for net upper and lower bounds
  BBox evalNetLowerBoundBB(unsigned netToEval) const;
  BBox evalNetUpperBoundBB(unsigned netToEval) const;
  void updateNetUpperAndLowerBounds(void);
  BBox getNetLowerBoundBB(unsigned netToEval) const {
    abkfatal(netToEval < _netLowerBounds.size(),
             "Attempt to get upper bound out of bounds");
    return _netLowerBounds[netToEval];
  }
  BBox getNetUpperBoundBB(unsigned netToEval) const {
    abkfatal(netToEval < _netUpperBounds.size(),
             "Attempt to get upper bound out of bounds");
    return _netUpperBounds[netToEval];
  }
  double getWorstCaseHPWL(void);
  double getBestCaseHPWL(void);
  // End functions added by DP for net upper and lower bounds

  void getBinMembership(std::vector<int>& nodeBin) const;
  // use the given bins to generate membership
  void getBinMembership(std::vector<int>& nodeBin,
                        const std::vector<CapoBin*>& bins) const;
  void printBinMembership(std::ostream& os) const;
  void printBinMembership(const char* fileName) const;

  void getHierCellNames(std::vector<std::string>& cellNames) const;
  void printHierCellNames(std::ostream& os) const;
  void printHierCellNames(const char* fileName) const;

  void saveSmallProblem(const PartitioningProblemForCapo& problem,
                        const CapoBin& bin);

  friend std::ostream& operator<<(std::ostream& os,
                                  std::vector<CapoBin*>& layout);

  // added by sadya to save the bins info as a floorplannning instance
  void saveBinsAsFloorplan(const char* fileName, std::vector<CapoBin*> bins);
  // save current snapshot of bins in a copy
  void saveBinsCopy();
  void solveQuadraticMinSoln();

  // get number of cells with region constraints in a bin
  unsigned getNumRegionedCellsInBin(CapoBin* bin);

  void boost(std::vector<BBox>& bboxes);
  void boostNets(CapoBin& bin, std::vector<BBox>&);
  void deboostNets(CapoBin&);
  double squareCost();
  void printNetStat();
  int probeBin(CapoBin& bin);
  void updateCell2Bin(const CapoBin* bin);
  void undoLayer(std::vector<CapoBin*>&, std::vector<CapoBin*>&);

  //  protected:
  // added by royj: necessary for integrating routing into the placement flow
  void capoPreMainLoop(unsigned& maxBins, unsigned& saveAtBins,
                       bool& haveUnplacedBins, bool& doVCuts, bool& saveBins);
  bool capoMainLoopTest(bool haveUnplacedBins, unsigned maxBins);
  void capoMainLoop(unsigned& saveAtBins, bool& haveUnplacedBins,
                    bool& saveBins, bool& doVCuts,
                    std::vector<CapoBin*>& newLayout,
                    std::vector<CapoBin*>& extraLayout,
                    std::vector<CapoBin*>& bestLayout);
  void capoPostMainLoop(Timer& capoTimer, RBPlacement& rbplace, bool& saveBins);

  void recordCutLine(const Point& a, const Point& b, unsigned bin1,
                     unsigned bin2);
  void recordCutLayer(unsigned layer);

  class FastPlaceMove {
   public:
    unsigned cellId;
    double cost;
    CapoBin *source, *dest;
    Point pos;

    FastPlaceMove()
        : cellId(UINT_MAX),
          cost(DBL_MAX),
          source(NULL),
          dest(NULL),
          pos(0., 0.) {}

    FastPlaceMove(const FastPlaceMove& m)
        : cellId(m.cellId),
          cost(m.cost),
          source(m.source),
          dest(m.dest),
          pos(m.pos) {}
  };

  class DecreasingByCost {
    const std::vector<FastPlaceMove>& _moveList;

   public:
    DecreasingByCost(const std::vector<FastPlaceMove>& moveList)
        : _moveList(moveList) {}

    bool operator()(unsigned a, unsigned b) {
      return (_moveList[a].cost > _moveList[b].cost);
    }
  };

  unsigned fastPlaceMovesMade;
  bool doFastPlaceMoves();
  fastSteiner::Point* _fastStPt;
  long* _fastStParent;
  double calculateSteinerFromPointSet(const std::vector<Point>& pointSet,
#ifdef USEFLUTE
                                      bool useFastSt, bool useFLUTE,
                                      bool useMST) const;
#else
                                      bool useFastSt, bool useMST) const;
#endif
  double calcFastPlaceCost(unsigned cellId, const Point& newCellPos,
                           std::vector<double>& curNetLens);
  void chooseFastPlaceMove(FastPlaceMove& move,
                           std::vector<double>& curNetLens);
  void makeFastPlaceMove(FastPlaceMove& move);
  void addCellsNeedingUpdate(CapoBin* bin, BitBoard& needsUpdate);
  void addCellsNeedingUpdate(unsigned cellId, BitBoard& needsUpdate);
  void invalidateNetLens(unsigned cellId, std::vector<double>& curNetLens);
#ifdef USEHMETIS
  void turnOffHMetis(void) { _params.splitterParams.useHMetis = false; }
#endif

  void randomlyMoveCellsInTheirBins(void);
};

#endif
