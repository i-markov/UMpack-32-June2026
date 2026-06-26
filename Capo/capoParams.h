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

// Created: Feb 15, 2000 by Andrew Caldwell from capoPlacer.h

#ifndef __CAPOPARAMS_H__
#define __CAPOPARAMS_H__

#include <iosfwd>

#include "ABKCommon/paramproc.h"
#include "ABKCommon/infolines.h"
#include <string>
#include <vector>
#include "ABKCommon/verbosity.h"
#include "MLPart/mlPart.h"
#include "ShellPart/bbFMPart.h"
#include "SmallPlacers/baseSmallPl.h"

class CapoSplitterParams {
 public:
  // if doRepartitioning, then the first tolerance is always
  // 20%, and the second tolerance follows the parameters below
  bool doRepartitioning;

  // if repartSmallWS then whenever partitioning a bin with < 2% WhiteSpace
  // use Repartitioning
  bool repartSmallWS;

  // if repartShiftCL then whenever cut-line shifts we repart, with smaller
  // tolerance
  bool repartShiftCL;

  // if useQuadCluster then use quadratic Clustering during MLPart
  bool useQuadCluster;

  // if <useWSTolMethod>, then use the computation from the
  // whitespace report for horizontal cuts.  Otherwise,
  // use <constantTolerance>.
  // If doRepartitioning is enabled, then v-cuts will
  // repartition with 1% tolerance after shifting.
  bool useWSTolMethod;
  double defaultTolerance;
  double defaultRepartTolerance;
  double minTolerance;
  double minRepartTolerance;
  double tolCap;

  // now both H-cuts and V-cuts use useWSTolMethod.
  // if uniformWS is enabled then V-cuts always use a tolerance of
  //<constantTolerance>, and shift the cutline.

  double safeWS;
  double minLocalWS;
  double smallWS;

  bool congestionShifting;
  bool alwaysCongShift;
  double congExponent;

  bool WSOracle;
  double WSOracleTol;

  bool eco;
  bool ecoShift;
#ifdef USEFLOW
  bool flowECO;
#endif
  double ecoOverfillPct;
  double ecoImprovePct;

  unsigned numMLSets;

  bool useNameHierarchy;
  std::string delimiters;

  // by sadya in SYNP to honour region constraints
  bool useRgnConstr;

  // <aaronnn> redistribute WS after partitioning by utilization density
  bool adjustWSByDensity;

  // by sadya in SYNP to try different partition capacities in presence of
  // filler cells to force the std-cells in the center of the design
  bool fillerCellsPartCaps;

  Verbosity verb;

  // by royj for adjusting branching factor
  double bFactorAdjustCutoff;

  // by royj for adaptive branching
  double extraBranchCutoff;

  // by royj for saving memory on partitionings
  double termMergeThreshold;

  bool newFuzzy;
  bool newFuzzyOnlyRepart;

  // by DAP for cutline planning
  bool useCutOracle;
  bool useOracleDirection;
  bool useACG;
  bool tightenACG;
#ifdef USEHMETIS
  bool useHMetis;
  double hmetisProblemThresh;
#endif
  bool useFMPartPlus;

  unsigned numCellGroups;
  unsigned cellGroupsCutoff;

  bool netDegreeWeighting;
  double netDegreeWeightExponent;

  CapoSplitterParams(int argc, const char *argv[]);
  CapoSplitterParams(Verbosity verbosity = Verbosity("0_0_0"))
      : doRepartitioning(false),
        repartSmallWS(true),
        repartShiftCL(true),
        useWSTolMethod(true),
        defaultTolerance(0.20),
        defaultRepartTolerance(0.10),
        minTolerance(0.05),
        minRepartTolerance(0.02),
        tolCap(1.),
        safeWS(0.05),
        minLocalWS(0.025),
        smallWS(0.02),
        congestionShifting(false),
        alwaysCongShift(true),
        congExponent(1.),
        WSOracle(false),
        WSOracleTol(0.02),
        eco(false),
        ecoShift(false),
#ifdef USEFLOW
        flowECO(false),
#endif
        ecoOverfillPct(0.),
        ecoImprovePct(1.),
        numMLSets(2),
        useNameHierarchy(false),
        delimiters("/"),
        useRgnConstr(true),
        adjustWSByDensity(false),
        fillerCellsPartCaps(false),
        verb(verbosity),
        bFactorAdjustCutoff(0.),
        extraBranchCutoff(4.),
        termMergeThreshold(0.),
        newFuzzy(false),
        newFuzzyOnlyRepart(false),
        useCutOracle(false),
        useOracleDirection(false),
        useACG(false),
        tightenACG(false),
#ifdef USEHMETIS
        useHMetis(true),
        hmetisProblemThresh(0.75),
#endif
        useFMPartPlus(false),
        numCellGroups(1),
        cellGroupsCutoff(35),
        netDegreeWeighting(false),
        netDegreeWeightExponent(0.) {
  }

  CapoSplitterParams(const CapoSplitterParams &srcParams);

  friend std::ostream &operator<<(std::ostream &,
                                  const CapoSplitterParams &params);
};

enum SeedPlacerType { DumbPlacerSeed, BaryPlacerSeed, WeiszPlacerSeed };
enum ReplaceSmallBinsType { Never, AtTheEnd, AtEveryLayer };

class CapoParameters {
  void setDefaults();

 public:
  Verbosity verb;

  MaxMem *maxMem;

  //  PARAMETERS FOR TOP-DOWN FLOW
  unsigned stopAtBins;  // if !=0, stop when >= this many bins
  // default == 0 => run till end.

  unsigned saveAtBins;  // if !=0, save when >= this many bins
  // default == 0 => run till end.

  ReplaceSmallBinsType replaceSmallBins;

  bool useActualPinLocs;

  // from royj: if automatic parameter is false, then Capo initialzes
  // variables and exits before doing anything so that another mechanism
  // can decide how to use the Capo utilities (i.e. necessary for
  // simultaneous placement and routing)
  bool automatic;

  //   PARAMETERS EQUALLY APPLICABLE TO EACH LEVEL
  unsigned smallPartThreshold;  // for problems of size <= the
  // threshold, use smallPartitioner
  unsigned smallPlaceThreshold;  // call the end-case placer on
  // one-row problems with <= threshold #cells
  unsigned smallSplitThreshold;  // call the small bin splitter on
                                 // bins with <= threshold # cells

  bool doOverlapping;
  unsigned startOverlappingLayer;  // the first level to overlap on
  unsigned endOverlappingLayer;    // the last level to overlap at

  // parameters for capo boosting. by sherief. ported by sadya
  unsigned boostLayer;
  unsigned boostFactor;
  bool boost;
  bool feedback;
  bool acc;  // accelerated feedback

  // parameters for feedback by Sherief.
  unsigned termIter;
  unsigned contMode;

  // by DAP to turn on Weighted Terminal Propagation based partitioning
  bool wtCut;
  // turn on Weighted Terminal Prop on bins larger than a certain size
  unsigned wtCutThreshold;
  // by DAP to trim bboxes and set cells by default at the center of gravity
  bool cogLocation;

  unsigned samples;

  // <aaronnn>: ThunderPart (global placement refinement - optimizes for HPWL)
  bool thunder;
  double thunderMoveMultiplier;

  // <aaronnn>: SCAMPI activation threshold (num macros)
  unsigned scampiThresh;

  // <aaronnn>: enable SOR for floorplanning certain instances
  bool useAnalytFP;

  // jflu: use cell bloating for whitespace
  bool useBloating;
  double bloatCongPercentile;

  // jflu: divide the net weights by square of bin utilization
  bool WSadjustment;
  double termWeightModifier;

  // jflu: turn on steiner tree evaluation methods
  bool useHPWL;
  bool useFastSt;
#ifdef USEFLUTE
  bool useFLUTE;
#endif
  bool useMST;

  bool doFastPlaceMoves;
  bool fastPlaceHPWL;
  bool fastPlaceFastSt;
#ifdef USEFLUTE
  bool fastPlaceFLUTE;
#endif

  // <aaronnn>: param for layer seeding
  std::vector<unsigned> layerSeeds;

  // jflu & aaronnn: use intermediate placement to fix cells before
  // partitioning?
  bool prePartFix;
  double prePartFixMargin;
  double prePartFixLimit;

  bool prePartSavePl;

  bool selectiveFloorplanning;

  // <aaronnn>: macro clustering
  bool clusterMacros;

  // <aaronnn>: FP obstacle handling
  bool FPEvadeObstacles;
  // sadya: FP obstacle handling threshold
  double FPEvadeObstaclesThresh;

  // <aaronnn> try FPing after partitioning. If fails, go up and try again
  bool lookAheadFP;

  // by sadya in SYNP to honour group region constraints
  bool useGrpConstr;

  // by sadya in SYNP to use predetermined cutlines while performing placement
  // Input is from a file called preDeterminedCultines.txt in the run directory
  bool usePredeterminedCutlines;

  // by sadya in SYNP to use determine cutlines aligned to big obstacles edges
  bool alignCutline2Obs;

  //   PARTITIONER PARAMETERS
  MLPartParams mlParams;
  SmallPlParams smplParams;
  BBFMPart::Parameters bbfmParams;

  // by sadya for top down congestion driven partitioning
  bool tdCongestionCtl;
  // by sadya to incorporate net effect
  bool capoNE;
  bool noCOG;    // don't use center of gravity constraints
  bool useQuad;  // generate quadratic min soln before each layer. usefull for
                 // terminal propagation
  bool useLinearQP;  // use gordian-l style linearization during QP solutions
  unsigned
      quadSkipNetsLargerThan;  // during SOR for quadratic min soln, skip nets
                               // larger than this number for time efficiency

  bool mixedSize;
  bool plotFPOutlines;
  bool plotCutLines;
  unsigned cutLineLayerCap;

  std::string FPOutlinesFileName;
  std::string CutLinesFileName;

  std::string netStatsFileName;

  // lookAhead > 1 requires use of the LASplitter
  unsigned lookAhead;

  bool allowRowSplits;
  CapoSplitterParams splitterParams;

  // for mixed-size placement

  //  INFO-TRACKING PARAMETERS
  bool plotBins;
  bool saveLayerBBoxes;
  bool saveBins;
  bool saveSmallPlProb;

  bool saveBinsFloorplan;  // save bins before detailed placement

  std::string FPrep;  // by royj for switching the floorplanning
                       // representation, currently either "SeqPair" or "BTree"
  bool noParquetRotation;
  bool snapMacrosX;
  bool snapMacrosY;
  bool removeMacroOverlap;
  bool allMacroSoft;  // by hhchan to mark all macros as soft blocks
  double maxAR;       // by hhchan, max aspect ratio for soft blocks
  double minAR;       // '' min aspect ratio ''
#ifdef USEPATOMA
  bool usePatoma;  // hhchan: to use Patoma instead of standard FP-flow
#endif
#ifdef USEFLOW
  bool flowAfterParquet;
  double flowLimitRatio;
#endif
  std::string softFileName;  // hhchan: determine which blocks are soft/hard

#ifdef USETRAFFIC
  bool useTraffic;
  double trafficOdds, trafficThresh;
#endif

  double vertFactor;
  double ARStretchFactor;  // by sadya in SYNP. determing X/Y HPWL ratio by
                           // altering cut direction based on AR of bin

  bool annealToCutLines;
  double parquetShrink;

  double keepHGraphFraction;  // keep the HGraph in memory during ML
                              // partitioning if the new HGraph will have fewer
                              // than this fraction of the number of nodes

  CapoParameters();
  CapoParameters(int argc, const char *argv[], bool saveMem);
  CapoParameters(const CapoParameters &srcParams);
  ~CapoParameters();

  void printHelp();

  friend std::ostream &operator<<(std::ostream &os, CapoParameters &params);
};

#endif
