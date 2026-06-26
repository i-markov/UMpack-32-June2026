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

// Created: 98/02/15 by Andrew Caldwell (caldwell@cs.ucla.edu)

// CHANGES
// 980215 aec created capoParams.cxx from capoPlacer.cxx
// 980313 ilm added CapoPlacer::Parameters ctor from argc, argv
//            reworked CapoPlacer::Parameters and its op<<
// 980325 ilm added n,m and partFuzziness to CapoPlacer::Parameters
// 980401 ilm added savePartProb parameter

#ifdef _MSC_VER
#pragma warning(disable : 4786)
#endif

#include "capoParams.h"

#include <sstream>

#include "ABKCommon/abkassert.h"
#include "HGraphWDims/hgWDims.h"

using std::cout;
using std::endl;
using std::istringstream;
using std::max;
using std::min;
using std::ostream;
using std::string;

void CapoParameters::setDefaults() {
  verb = Verbosity("1 1 1");
  maxMem = NULL;

  stopAtBins = 0;
  saveAtBins = 0;
  mlParams.netThreshold = 0;  // don't threshold

  plotBins = false;
  boost = false;
  saveLayerBBoxes = false;
  saveBins = false;
  saveSmallPlProb = false;
  saveBinsFloorplan = false;

  smallPartThreshold = 35;
  smallPlaceThreshold = 7;
  smallSplitThreshold = 200;

  boostLayer = 6;
  boostFactor = 2;
  termIter = 3;
  contMode = 1;
  feedback = true;
  acc = true;

  wtCut = true;
  wtCutThreshold = 0;
  cogLocation = false;
  samples = 0;

  // <aaronnn>
  thunder = false;
  thunderMoveMultiplier = 1;
  layerSeeds.clear();

  // <aaronnn>
  scampiThresh = 100;

  // <aaronnn>
  useAnalytFP = false;

  // jflu
  useBloating = false;
  bloatCongPercentile = 0.95;
  WSadjustment = true;
  termWeightModifier = 2;

  useHPWL = true;
  useFastSt = false;
#ifdef USEFLUTE
  useFLUTE = false;
#endif
  useMST = false;

  doFastPlaceMoves = false;
  fastPlaceHPWL = true;
  fastPlaceFastSt = false;
#ifdef USEFLUTE
  fastPlaceFLUTE = false;
#endif

  // jflu & aaronnn
  prePartFix = false;
  prePartFixMargin = 0.2;
  prePartFixLimit = 0.1;
  prePartSavePl = false;

  // only certain macros will be placed
  // after floorplanning
  selectiveFloorplanning = false;

  // <aaronnn> clusterMacros
  clusterMacros = true;

  // <aaronnn> obstacle handling in FP
  FPEvadeObstacles = true;
  FPEvadeObstaclesThresh = 0;

  // <aaronnn> perform lookahead FP to detect FP failures and go up if necessary
  lookAheadFP = true;

  doOverlapping = false;
  startOverlappingLayer = 0;
  endOverlappingLayer = UINT_MAX;

  replaceSmallBins = AtTheEnd;
  //   replaceSmallBins = Never;
  useActualPinLocs = true;

  tdCongestionCtl = false;
  capoNE = false;
  noCOG = false;
  useQuad = false;
  useLinearQP = false;
  quadSkipNetsLargerThan = 150;
  mixedSize = true;
  plotFPOutlines = false;
  plotCutLines = false;
  cutLineLayerCap = UINT_MAX;

  FPOutlinesFileName = "floorplannedBins.dat";
  CutLinesFileName = "cutLines.dat";
  netStatsFileName = "";

  lookAhead = 0;
  allowRowSplits = false;
  useGrpConstr = true;
  usePredeterminedCutlines = true;
  alignCutline2Obs = false;

  // added by royj
  automatic = true;

  FPrep = "SeqPair";
  noParquetRotation = false;
  snapMacrosX = false;
  snapMacrosY = false;
  removeMacroOverlap = true;
  allMacroSoft = false;
  maxAR = HGraphWDimensions::DEFAULT_MAX_AR;
  minAR = HGraphWDimensions::DEFAULT_MIN_AR;
#ifdef USEPATOMA
  usePatoma = false;
#endif
#ifdef USEFLOW
  flowAfterParquet = false;
  flowLimitRatio = 0.01;
#endif
  softFileName = "";

  vertFactor = 7.0;
  ARStretchFactor = -1.;

  annealToCutLines = false;
  parquetShrink = -1.;

  keepHGraphFraction = 2.;

#ifdef USETRAFFIC
  useTraffic = false;
  trafficOdds = 50.;
  trafficThresh = .1;
#endif
}

CapoParameters::CapoParameters() : verb("1 1 1") { setDefaults(); }

CapoParameters::CapoParameters(const CapoParameters &srcParams)
    : verb(srcParams.verb),
      maxMem(srcParams.maxMem),
      stopAtBins(srcParams.stopAtBins),
      saveAtBins(srcParams.saveAtBins),
      replaceSmallBins(srcParams.replaceSmallBins),
      useActualPinLocs(srcParams.useActualPinLocs),
      automatic(srcParams.automatic),
      smallPartThreshold(srcParams.smallPartThreshold),
      smallPlaceThreshold(srcParams.smallPlaceThreshold),
      smallSplitThreshold(srcParams.smallSplitThreshold),
      doOverlapping(srcParams.doOverlapping),
      startOverlappingLayer(srcParams.startOverlappingLayer),
      endOverlappingLayer(srcParams.endOverlappingLayer),
      boostLayer(srcParams.boostLayer),
      boostFactor(srcParams.boostFactor),
      boost(srcParams.boost),
      feedback(srcParams.feedback),
      acc(srcParams.acc),
      termIter(srcParams.termIter),
      contMode(srcParams.contMode),
      wtCut(srcParams.wtCut),
      wtCutThreshold(srcParams.wtCutThreshold),
      cogLocation(srcParams.cogLocation),
      samples(srcParams.samples),
      thunder(srcParams.thunder),
      thunderMoveMultiplier(srcParams.thunderMoveMultiplier),
      scampiThresh(srcParams.scampiThresh),
      useAnalytFP(srcParams.useAnalytFP),
      useBloating(srcParams.useBloating),
      bloatCongPercentile(srcParams.bloatCongPercentile),
      WSadjustment(srcParams.WSadjustment),
      termWeightModifier(srcParams.termWeightModifier),
      useHPWL(srcParams.useHPWL),
      useFastSt(srcParams.useFastSt),
#ifdef USEFLUTE
      useFLUTE(srcParams.useFLUTE),
#endif
      useMST(srcParams.useMST),
      doFastPlaceMoves(srcParams.doFastPlaceMoves),
      fastPlaceHPWL(srcParams.fastPlaceHPWL),
      fastPlaceFastSt(srcParams.fastPlaceFastSt),
#ifdef USEFLUTE
      fastPlaceFLUTE(srcParams.fastPlaceFLUTE),
#endif
      layerSeeds(srcParams.layerSeeds),
      prePartFix(srcParams.prePartFix),
      prePartFixMargin(srcParams.prePartFixMargin),
      prePartFixLimit(srcParams.prePartFixLimit),
      prePartSavePl(srcParams.prePartSavePl),
      selectiveFloorplanning(srcParams.selectiveFloorplanning),
      clusterMacros(srcParams.clusterMacros),
      FPEvadeObstacles(srcParams.FPEvadeObstacles),
      FPEvadeObstaclesThresh(srcParams.FPEvadeObstaclesThresh),
      lookAheadFP(srcParams.lookAheadFP),
      useGrpConstr(srcParams.useGrpConstr),
      usePredeterminedCutlines(srcParams.usePredeterminedCutlines),
      alignCutline2Obs(srcParams.alignCutline2Obs),
      mlParams(srcParams.mlParams),
      smplParams(srcParams.smplParams),
      bbfmParams(srcParams.bbfmParams),
      tdCongestionCtl(srcParams.tdCongestionCtl),
      capoNE(srcParams.capoNE),
      noCOG(srcParams.noCOG),
      useQuad(srcParams.useQuad),
      useLinearQP(srcParams.useLinearQP),
      quadSkipNetsLargerThan(srcParams.quadSkipNetsLargerThan),
      mixedSize(srcParams.mixedSize),
      plotFPOutlines(srcParams.plotFPOutlines),
      plotCutLines(srcParams.plotCutLines),
      cutLineLayerCap(srcParams.cutLineLayerCap),
      FPOutlinesFileName(srcParams.FPOutlinesFileName),
      CutLinesFileName(srcParams.CutLinesFileName),
      netStatsFileName(srcParams.netStatsFileName),
      lookAhead(srcParams.lookAhead),
      allowRowSplits(srcParams.allowRowSplits),
      splitterParams(srcParams.splitterParams),
      plotBins(srcParams.plotBins),
      saveLayerBBoxes(srcParams.saveLayerBBoxes),
      saveBins(srcParams.saveBins),
      saveSmallPlProb(srcParams.saveSmallPlProb),
      saveBinsFloorplan(srcParams.saveBinsFloorplan),
      FPrep(srcParams.FPrep),
      noParquetRotation(srcParams.noParquetRotation),
      snapMacrosX(srcParams.snapMacrosX),
      snapMacrosY(srcParams.snapMacrosY),
      removeMacroOverlap(srcParams.removeMacroOverlap),
      allMacroSoft(srcParams.allMacroSoft),
      maxAR(srcParams.maxAR),
      minAR(srcParams.minAR),
#ifdef USEPATOMA
      usePatoma(srcParams.usePatoma),
#endif
#ifdef USEFLOW
      flowAfterParquet(srcParams.flowAfterParquet),
      flowLimitRatio(srcParams.flowLimitRatio),
#endif
      softFileName(srcParams.softFileName),
#ifdef USETRAFFIC
      useTraffic(srcParams.useTraffic),
      trafficOdds(srcParams.trafficOdds),
      trafficThresh(srcParams.trafficThresh),
#endif
      vertFactor(srcParams.vertFactor),
      ARStretchFactor(srcParams.ARStretchFactor),
      annealToCutLines(srcParams.annealToCutLines),
      parquetShrink(srcParams.parquetShrink),
      keepHGraphFraction(srcParams.keepHGraphFraction) {
}

CapoParameters::CapoParameters(int argc, const char *argv[], bool saveMem)
    : mlParams(argc, argv),
      smplParams(argc, argv),
      bbfmParams(argc, argv),
      splitterParams(argc, argv) {
  setDefaults();

  mlParams.saveMem = saveMem;
  bbfmParams.fmParams.saveMem = saveMem;

  BoolParam helpRequest0("help", argc, argv);
  BoolParam helpRequest1("h", argc, argv);

  if (helpRequest0.found() || helpRequest1.found()) {
    printHelp();
    return;
  }

  //----------------------------------------------------------------

  StringParam verb_("verb", argc, argv);

  // general Capo params..
  UnsignedParam stopAtBins_("stopAtBins", argc, argv);
  UnsignedParam saveAtBins_("saveAtBins", argc, argv);
  UnsignedParam smallPartThreshold_("smallPartThreshold", argc, argv);
  UnsignedParam smallPlaceThreshold_("smallPlaceThreshold", argc, argv);
  UnsignedParam smallSplitThreshold_("smallSplitThreshold", argc, argv);
  UnsignedParam boostLayer_("boostLayer", argc, argv);
  UnsignedParam boostFactor_("boostFactor", argc, argv);
  UnsignedParam termIter_("termIter", argc, argv);
  UnsignedParam contMode_("contMode", argc, argv);
  BoolParam noWtCut("noWtCut", argc, argv);
  UnsignedParam wtCutThresh("wtCutThresh", argc, argv);
  BoolParam cogLocation_("cogLocation", argc, argv);
  UnsignedParam samples_("samples", argc, argv);
  StringParam replaceSmallBins_("replaceSmallBins", argc, argv);
  BoolParam useActualPinLocs_("useActualPinLocs", argc, argv);
  BoolParam tdCongestionCtl_("tdCongestionCtl", argc, argv);
  BoolParam capoNE_("capoNE", argc, argv);
  BoolParam noCOG_("noCOG", argc, argv);
  BoolParam useQuad_("useQuad", argc, argv);
  BoolParam useLinearQP_("useLinearQP", argc, argv);
  UnsignedParam quadSkipNetsLargerThan_("quadSkipNetsLargerThan", argc, argv);
  BoolParam noMS_("noMS", argc, argv);
  BoolParam useSeqPair_("useSeqPair", argc, argv);
  BoolParam useBstarTree_("useBstarTree", argc, argv);
  BoolParam allMacroSoft_("allMacroSoft", argc, argv);
  DoubleParam maxAR_("maxAR", argc, argv);
  DoubleParam minAR_("minAR", argc, argv);
#ifdef USEPATOMA
  BoolParam usePatoma_("usePatoma", argc, argv);
#endif
#ifdef USEFLOW
  BoolParam flowAfterParquet_("flowAfterParquet", argc, argv);
  DoubleParam flowLimitRatio_("flowLimitPct", argc, argv);
#endif
  StringParam softFileName_("selectSoft", argc, argv);
  DoubleParam vertFactor_("vertFactor", argc, argv);
  DoubleParam ARStretchFactor_("ARStretchFactor", argc, argv);
  BoolParam annealToCutLines_("annealToCutLines", argc, argv);
  DoubleParam parquetShrink_("parquetShrink", argc, argv);

  // Partitioner params...
  UnsignedParam lookAhead_("lookAhead", argc, argv);
  BoolParam allowRowSplits_("allowRowSplits", argc, argv);

  // <aaronnn> ThunderPart params...
  BoolParam thunder_("thunder", argc, argv);
  BoolParam noThunder_("noThunder", argc, argv);
  DoubleParam thunderMoveMultiplier_("thunderMoveMultiplier", argc, argv);

  // <aaronnn>
  UnsignedParam scampiThresh_("scampiThresh", argc, argv);

  // <aaronnn>
  BoolParam useAnalytFP_("useAnalytFP", argc, argv);

  // jflu: use to plot congestion estimation based on connectivity at each major
  // layer
  BoolParam useBloating_("useBloating", argc, argv);
  DoubleParam bloatCongPercentile_("bloatCongPercentile", argc, argv);

  // jflu: use some power of utilization to weight terminal bound nets
  BoolParam noWSadjustment("noWSadjustment", argc, argv);
  DoubleParam useTermWeightModifier("termWeightModifier", argc, argv);

  BoolParam useHPWL_("useHPWL", argc, argv);
  BoolParam useFastSt_("useFastSt", argc, argv);
#ifdef USEFLUTE
  BoolParam useFLUTE_("useFLUTE", argc, argv);
#endif
  BoolParam useMST_("useMST", argc, argv);

  BoolParam doFastPlaceMoves_("doFastPlaceMoves", argc, argv);
  BoolParam fastPlaceHPWL_("fastPlaceHPWL", argc, argv);
  BoolParam fastPlaceFastSt_("fastPlaceFastSt", argc, argv);
#ifdef USEFLUTE
  BoolParam fastPlaceFLUTE_("fastPlaceFLUTE", argc, argv);
#endif

  // jflu & aaronnn pre-partitioning cell fixing params
  BoolParam prePartFix_("prePartFix", argc, argv);
  BoolParam noPrePartFix_("noPrePartFix", argc, argv);
  DoubleParam prePartFixMargin_("prePartFixMargin", argc, argv);
  DoubleParam prePartFixLimit_("prePartFixLimit", argc, argv);

  // aaronnn save intermediate placement before partitioning
  BoolParam prePartSavePl_("prePartSavePl", argc, argv);

  // <aaronnn> layer seed params
  StringParam setLayerSeeds_("setLayerSeeds", argc, argv);

  BoolParam selectiveFloorplanning_("selectiveFloorplanning", argc, argv);

  // <aaronnn> macro clustering params
  BoolParam noClusterMacros_("noClusterMacros", argc, argv);

  // <aaronnn> obstacle evasion params
  BoolParam noFPEvadeObstacles_("noFPEvadeObstacles", argc, argv);
  BoolParam FPEvadeObstaclesThresh_("FPEvadeObstaclesThresh", argc, argv);

  // <aaronnn> lookahead FP params
  BoolParam lookAheadFP_("lookAheadFP", argc, argv);

  // Optional info-tracking params
  BoolParam plotBins_("plotBins", argc, argv);
  BoolParam boost_("boost", argc, argv);
  BoolParam feedback_("feedback", argc, argv);
  BoolParam noFeedback_("noFeedback", argc, argv);
  // BoolParam	acc_       ("acc", argc, argv);
  BoolParam saveLayerBBoxes_("saveLayerBBoxes", argc, argv);
  BoolParam saveBins_("saveBins", argc, argv);
  BoolParam saveSmallPlProb_("saveSmallPlProb", argc, argv);
  BoolParam saveBinsFloorplan_("saveBinsFloorplan", argc, argv);
  BoolParam plotFPOutlines_("plotFPOutlines", argc, argv);
  BoolParam plotCutLines_("plotCutLines", argc, argv);
  UnsignedParam cutLineLayerCap_("cutLineLayerCap", argc, argv);

  StringParam logNetStats("logNetStats", argc, argv);

  // Overlapping Parameters
  BoolParam doOverlapping_("doOverlapping", argc, argv);
  UnsignedParam startOverlappingLayer_("startOverlappingLayer", argc, argv);
  UnsignedParam endOverlappingLayer_("endOverlappingLayer", argc, argv);

  // To handle group constraints in global placement
  BoolParam noGrpConstr_("noGrpConstr", argc, argv);
  BoolParam noPDC_("noPredeterminedCutlines", argc, argv);
  BoolParam alignCutline2Obs_("alignCutline2Obs", argc, argv);

  // meta-options
  BoolParam ispd05_("ispd05", argc, argv);
  BoolParam ispd06_("ispd06", argc, argv);
  BoolParam tryharder_("tryHarder", argc, argv);
  BoolParam faster_("faster", argc, argv);
  BoolParam routable_("routable", argc, argv);
  BoolParam ROOSTER("ROOSTER", argc, argv);
  BoolParam FLROOSTER("FLROOSTER", argc, argv);
  BoolParam SCAMPI("SCAMPI", argc, argv);

  //----------------------------------------------------------------

  // process meta-options first so they can be tweaked with other options
  if (routable_.found() || ROOSTER.found()) {
    wtCut = true;
    useHPWL = false;
    useFastSt = true;
    acc = false;
    feedback = true;
    if (routable_.found()) {
      cout << "Meta-Option \"routable\" detected:" << endl;
    }
    if (ROOSTER.found()) {
      cout << "Meta-Option \"ROOSTER\" detected:" << endl;
    }
    cout << "  - full feedback enabled" << endl;
    cout << "  - FastSt Steiner weighting enabled" << endl << endl;
  }

#ifdef USEFLUTE
  if (FLROOSTER.found()) {
    wtCut = true;
    useHPWL = false;
    useFLUTE = true;
    acc = false;
    feedback = true;
    cout << "Meta-Option \"FLROOSTER\" detected:" << endl;
    cout << "  - full feedback enabled" << endl;
    cout << "  - FLUTE Steiner weighting enabled" << endl << endl;
  }
#endif

  if (ispd05_.found()) {
    // no rotations for parquet
    noParquetRotation = true;
    // snap macros to rows/sites
    snapMacrosX = true;
    snapMacrosY = true;
    cout << "Meta-Option \"ispd05\" detected:" << endl;
    cout << "  - macros will not be rotated, flipped, or mirrored" << endl;
    cout << "  - macros will be snapped to row and site lines" << endl << endl;
  }

  if (ispd06_.found()) {
    wtCut = false;
    useHPWL = false;
    // no rotations for parquet
    noParquetRotation = true;
    // snap macros to rows/sites
    snapMacrosX = true;
    snapMacrosY = true;
    // disable feedback
    acc = false;
    feedback = false;
    cout << "Meta-Option \"ispd06\" detected:" << endl;
    cout << "  - macros will not be rotated, flipped, or mirrored" << endl;
    cout << "  - macros will be snapped to row and site lines" << endl;
    cout << "  - weighted terminal propagation disabled" << endl;
    cout << "  - feedback disabled" << endl << endl;
  }

  if (faster_.found() || SCAMPI.found()) {
    wtCut = false;

    useHPWL = false;
    useFastSt = false;
#ifdef USEFLUTE
    useFLUTE = false;
#endif
    useMST = false;
    acc = false;
    feedback = false;
    if (faster_.found()) {
      cout << "Meta-Option \"faster\" detected:" << endl;
    }
    if (SCAMPI.found()) {
      cout << "Meta-Option \"SCAMPI\" detected:" << endl;
    }
    cout << "  - weighted cut disabled" << endl;
    cout << "  - feedback disabled" << endl;
  }

  if (tryharder_.found()) {
    acc = false;
    feedback = true;
    wtCut = true;
    useHPWL = true;
    cout << "Meta-Option \"tryHarder\" detected:" << endl;
    cout << "  - full feedback enabled" << endl << endl;
  }

  if (SCAMPI.found()) {
    FPrep = "BTree";
    cout << "  - using B*Tree as the floorplanning representation" << endl;
  }

  if (useFastSt_.found()) {
    wtCut = true;
    useFastSt = true;
    useHPWL = false;
  }

#ifdef USEFLUTE
  if (useFLUTE_.found()) {
    wtCut = true;
    useFLUTE = true;
    useHPWL = false;
  }
#endif

  if (useHPWL_.found()) {
    wtCut = true;
    useHPWL = true;
    useFastSt = false;
#ifdef USEFLUTE
    useFLUTE = false;
#endif
    useMST = false;
  }

  if (useMST_.found()) {
    wtCut = true;
    useHPWL = false;
    useFastSt = false;
#ifdef USEFLUTE
    useFLUTE = false;
#endif
    useMST = true;
  }

  if (noWtCut.found()) {
    wtCut = false;
    useHPWL = false;
    useFastSt = false;
#ifdef USEFLUTE
    useFLUTE = false;
#endif
    useMST = false;
  }

  if (samples_.found()) {
    samples = samples_;
  }

  if (smallPartThreshold_.found()) smallPartThreshold = smallPartThreshold_;
  if (smallPlaceThreshold_.found()) smallPlaceThreshold = smallPlaceThreshold_;
  if (smallSplitThreshold_.found()) smallSplitThreshold = smallSplitThreshold_;
  if (boostLayer_.found()) boostLayer = boostLayer_;
  if (boostFactor_.found()) boostFactor = boostFactor_;
  if (termIter_.found()) termIter = termIter_;
  if (contMode_.found()) contMode = contMode_;

  if (wtCutThresh.found()) {
    wtCutThreshold = wtCutThresh;
    wtCut = true;
  }

  if (cogLocation_.found()) cogLocation = true;

  if (replaceSmallBins_.found()) {
    if (!strcasecmp(replaceSmallBins_, "Never"))
      replaceSmallBins = Never;
    else if (!strcasecmp(replaceSmallBins_, "AtTheEnd"))
      replaceSmallBins = AtTheEnd;
    else if (!strcasecmp(replaceSmallBins_, "AtEveryLayer"))
      replaceSmallBins = AtEveryLayer;
    else
      abkfatal(0, "invalid string follows '-replaceSmallBins'");
  }

  if (useActualPinLocs_.found()) useActualPinLocs = true;

  if (verb_.found()) verb = Verbosity(verb_);

  if (stopAtBins_.found()) stopAtBins = stopAtBins_;
  if (saveAtBins_.found()) saveAtBins = saveAtBins_;

  if (noGrpConstr_.found()) useGrpConstr = false;
  if (noPDC_.found()) usePredeterminedCutlines = false;
  if (alignCutline2Obs_.found()) alignCutline2Obs = true;

  plotBins = plotBins_.found();
  boost = boost_.found();

  // feedback = feedback_.found();
  // acc = acc_.found();
  if (feedback_.found()) acc = false;

  if (noFeedback_.found()) {
    feedback = false;
    acc = false;
  }

  saveLayerBBoxes = saveLayerBBoxes_.found();
  saveBins = saveBins_.found();
  saveSmallPlProb = saveSmallPlProb_.found();
  saveBinsFloorplan = saveBinsFloorplan_.found();

  if (doOverlapping_.found()) doOverlapping = true;
  if (startOverlappingLayer_.found())
    startOverlappingLayer = startOverlappingLayer_;
  if (endOverlappingLayer_.found()) endOverlappingLayer = endOverlappingLayer_;

  if (tdCongestionCtl_.found()) tdCongestionCtl = tdCongestionCtl_;

  if (capoNE_.found()) capoNE = capoNE_;

  if (noCOG_.found()) noCOG = noCOG_;

  if (useQuad_.found()) useQuad = useQuad_;

  if (useLinearQP_.found()) useLinearQP = useLinearQP_;

  if (quadSkipNetsLargerThan_.found())
    quadSkipNetsLargerThan = quadSkipNetsLargerThan_;

  if (noMS_.found()) mixedSize = false;

  if (plotFPOutlines_.found()) plotFPOutlines = plotFPOutlines_;

  if (plotCutLines_.found()) plotCutLines = plotCutLines_;

  if (cutLineLayerCap_.found()) cutLineLayerCap = cutLineLayerCap_;

  if (logNetStats.found()) netStatsFileName = logNetStats;

  if (lookAhead_.found()) lookAhead = lookAhead_;

  if (allowRowSplits_.found()) allowRowSplits = true;

  if (useSeqPair_.found())
    FPrep = "SeqPair";
  else if (useBstarTree_.found())
    FPrep = "BTree";

  if (allMacroSoft_.found()) allMacroSoft = allMacroSoft_;

  if (maxAR_.found()) maxAR = maxAR_;
  if (minAR_.found()) minAR = minAR_;
#ifdef USEPATOMA
  if (usePatoma_.found()) usePatoma = true;
#endif
#ifdef USEFLOW
  if (flowAfterParquet_.found()) {
    flowAfterParquet = true;
    snapMacrosX = true;
    snapMacrosY = true;
  }

  if (flowLimitRatio_.found()) {
    flowAfterParquet = true;
    snapMacrosX = true;
    snapMacrosY = true;
    flowLimitRatio = 0.01 * flowLimitRatio_;
  }
#endif
  if (softFileName_.found()) softFileName = softFileName_;

  if (vertFactor_.found()) vertFactor = vertFactor_;

  if (ARStretchFactor_.found()) ARStretchFactor = ARStretchFactor_;

  if (annealToCutLines_.found()) annealToCutLines = annealToCutLines_;

  if (parquetShrink_.found()) {
    parquetShrink = parquetShrink_;
    if (greaterThanDouble(parquetShrink, 1.)) {
      parquetShrink *= 0.01;
    }
  }

  if (thunder_.found()) thunder = true;
  if (noThunder_.found()) thunder = false;
  if (thunderMoveMultiplier_.found())
    thunderMoveMultiplier = thunderMoveMultiplier_;

  if (scampiThresh_.found()) scampiThresh = scampiThresh_;

  if (useAnalytFP_.found()) useAnalytFP = useAnalytFP_;

  if (selectiveFloorplanning_.found()) selectiveFloorplanning = true;

  // <aaronnn> macro clustering
  if (noClusterMacros_.found()) clusterMacros = false;

  // <aaronnn> FP obstacle evasion
  if (noFPEvadeObstacles_.found()) FPEvadeObstacles = false;

  if (FPEvadeObstaclesThresh_.found())
    FPEvadeObstaclesThresh = 0.01 * FPEvadeObstaclesThresh_;

  // <aaronnn> FP after partitioning as a 'sanity check'
  if (lookAheadFP_.found()) lookAheadFP = true;

  if (useBloating_.found()) useBloating = true;

  if (bloatCongPercentile_.found())
    bloatCongPercentile = bloatCongPercentile_ / 100.;

  if (noWSadjustment.found()) WSadjustment = false;

  if (useTermWeightModifier.found()) termWeightModifier = useTermWeightModifier;

  if (setLayerSeeds_.found()) {
    // seeds were supplied for each layer - split ":" delimited string into
    // vector
    const char *setLayerSeeds = setLayerSeeds_;
    string toSplit(setLayerSeeds);
    string::size_type start = 0;
    string::size_type len = 0;
    string::size_type n;

    while ((n = toSplit.find(":", start)) != string::npos) {
      len = n - start;
      istringstream ss(toSplit.substr(start, len).c_str());
      unsigned seed;
      ss >> seed;
      // cout << "DBG val: " << seed << endl;
      layerSeeds.push_back(seed);
      start += len + 1;
    }
    // pick up the end (if there is one)
    if (start < toSplit.length()) {
      istringstream ss(toSplit.substr(start).c_str());
      unsigned seed;
      ss >> seed;
      // cout << "DBG val: " << seed << endl;
      layerSeeds.push_back(seed);
    }
  }

  if (prePartFix_.found()) prePartFix = true;
  if (noPrePartFix_.found()) prePartFix = false;
  if (prePartFixMargin_.found()) prePartFixMargin = prePartFixMargin_;
  if (prePartFixLimit_.found()) prePartFixLimit = prePartFixLimit_;

  if (prePartSavePl_.found()) prePartSavePl = prePartSavePl_;

  BoolParam useACG_("useACG", argc, argv);
  if (useACG_.found()) {
    useQuad = true;
    noCOG = true;
    cout << "Meta-Option \"useACG\" detected:" << endl;
    cout << "  - SOR turned on" << endl;
    cout << "  - SOR COG turned off" << endl;
  }

  BoolParam noCapoRemoveMacroOverlap("noCapoRemoveMacroOverlap", argc, argv);
  if (noCapoRemoveMacroOverlap.found()) {
    removeMacroOverlap = false;
  }

  if (doFastPlaceMoves_.found()) {
    doFastPlaceMoves = true;
  }

  if (fastPlaceHPWL_.found()) {
    doFastPlaceMoves = true;
    fastPlaceHPWL = true;
    fastPlaceFastSt = false;
#ifdef USEFLUTE
    fastPlaceFLUTE = false;
#endif
  }

  if (fastPlaceFastSt_.found()) {
    doFastPlaceMoves = true;
    fastPlaceHPWL = false;
    fastPlaceFastSt = true;
  }

#ifdef USEFLUTE
  if (fastPlaceFLUTE_.found()) {
    doFastPlaceMoves = true;
    fastPlaceHPWL = false;
    fastPlaceFLUTE = true;
  }
#endif

  UnsignedParam keepHGraphFraction_("keepHGraphFraction", argc, argv);

  if (keepHGraphFraction_.found()) {
    keepHGraphFraction = 0.01 * keepHGraphFraction_;
  }

#ifdef USETRAFFIC
  BoolParam useTraffic_("TRAFFIC", argc, argv);
  DoubleParam trafficOdds_("TRAFFICodds", argc, argv);
  DoubleParam trafficThresh_("TRAFFICthresh", argc, argv);

  if (useTraffic_.found()) {
    useTraffic = true;
    if (trafficOdds_.found()) trafficOdds = trafficOdds_;
    if (trafficThresh_.found()) trafficThresh = trafficThresh_ / 100.;
  }
#endif
}

ostream &operator<<(ostream &os, CapoParameters &params) {
  const char *tfStr[2] = {"false", "true"};

  os << "---- Capo Parameters     " << endl << endl;

  os << " StopAtBins:  ";
  if (params.stopAtBins)
    os << params.stopAtBins << endl;
  else
    os << " complete placement" << endl;
  os << " SaveAtBins:  ";
  if (params.saveAtBins)
    os << params.saveAtBins << endl;
  else
    os << " no bin save" << endl;

  os << " SmallPartThreshold:  " << params.smallPartThreshold << endl;
  os << " SmallPlaceThreshold: " << params.smallPlaceThreshold << endl;
  os << " SmallSplitThreshold: " << params.smallSplitThreshold << endl;

  os << " ReplaceSmallBins:  ";
  if (params.replaceSmallBins == Never)
    os << "Never" << endl;
  else if (params.replaceSmallBins == AtTheEnd)
    os << "AtTheEnd" << endl;
  else if (params.replaceSmallBins == AtEveryLayer)
    os << "AtEveryLayer" << endl;
  else
    abkfatal(0, "invalid value for replaceSmallBins");

  os << " Use actual pin locations: " << tfStr[params.useActualPinLocs] << endl;
  os << " Allow RowSplits :  " << tfStr[params.allowRowSplits] << endl;
  os << " Top-down congestion control :  " << tfStr[params.tdCongestionCtl]
     << endl;

  os << " Mixed-size placement : " << tfStr[params.mixedSize] << endl;
  os << " Selective floorplanning : " << tfStr[params.selectiveFloorplanning]
     << endl;
  os << " Use macro clustering : " << tfStr[params.clusterMacros] << endl;
  os << " Use obstacle evasion in floorplanning : "
     << tfStr[params.FPEvadeObstacles] << endl;
  if (params.FPEvadeObstacles)
    os << "   Obstacle evasion threshold: "
       << 100. * params.FPEvadeObstaclesThresh << "%" << endl;
  os << " Use look-ahead floorplanning : " << tfStr[params.lookAheadFP] << endl;
  os << " Plot Floorplan Outlines : " << tfStr[params.plotFPOutlines] << endl;
  if (params.plotFPOutlines) {
    os << "   Floorplan Outlines data filename : " << params.FPOutlinesFileName
       << endl;
  }
  os << " Plot Cut Lines : " << tfStr[params.plotCutLines] << endl;
  if (params.plotCutLines) {
    os << "   Layer cap for cut line plotting : " << params.cutLineLayerCap
       << endl;
    os << "   Cut Lines data filename : " << params.CutLinesFileName << endl;
  }

  if (params.netStatsFileName != "") {
    os << " Will log net statistics to file: \"" << params.netStatsFileName
       << "\"" << endl;
  }

  os << " Use boosting :         " << tfStr[params.boost] << endl;
  os << "     Boosting layer :         " << params.boostLayer << endl;
  os << "     Boosting factor :        " << params.boostFactor << endl;

  os << " Use feedback :         " << tfStr[params.feedback] << endl;
  os << "     termIter: " << params.termIter << endl;
  os << "     contMode: " << params.contMode << endl;
  os << "     accelerated FB: " << tfStr[params.acc] << endl;

  os << " Use ThunderPart :      " << tfStr[params.thunder] << endl;
  os << "     moveMultiplier : " << params.thunderMoveMultiplier << endl;

  os << " Use intermediate placements to fix cells before partitioning : "
     << tfStr[params.prePartFix] << endl;
  os << "     fixMargin : " << 100. * params.prePartFixMargin << "%" << endl;
  os << "     fixLimit  : " << 100. * params.prePartFixLimit << "%" << endl;

  bool usingSteiner = params.useFastSt;
#ifdef USEFLUTE
  usingSteiner = usingSteiner || params.useFLUTE;
#endif

  os << " Steiner Tree Evaluator Used: " << tfStr[usingSteiner] << endl;
  if (usingSteiner) {
    os << "     Fast Steiner: " << tfStr[params.useFastSt] << endl;
#ifdef USEFLUTE
    os << "     FLUTE       : " << tfStr[params.useFLUTE] << endl;
#endif
  }

  os << " Using MST for weighted terminal propagation: " << tfStr[params.useMST]
     << endl;

  os << " Do Fast Place Movement after regular Capo layers : "
     << tfStr[params.doFastPlaceMoves] << endl;
  if (params.doFastPlaceMoves) {
    os << "     Use HPWL        : " << tfStr[params.fastPlaceHPWL] << endl;
    os << "     Use Fast Steiner: " << tfStr[params.fastPlaceFastSt] << endl;
#ifdef USEFLUTE
    os << "     Use FLUTE       : " << tfStr[params.fastPlaceFLUTE] << endl;
#endif
  }

  if (!params.layerSeeds.empty()) {
    os << " Set layer seeds : ";
    for (unsigned i = 0; i < params.layerSeeds.size(); i++) {
      if (i != 0) os << ":";
      os << params.layerSeeds[i];
    }
    os << endl;
  }

  os << " Using Cell Bloating : " << tfStr[params.useBloating] << endl;
  if (params.useBloating)
    os << "  Bloating Congestion Percentile : "
       << params.bloatCongPercentile * 100. << "%" << endl;
  os << " Weighted Terminal Propagation :  " << tfStr[params.wtCut] << endl;
  if (params.wtCut) {
    os << "   Node Threshold for wtCut : " << params.wtCutThreshold << endl;
    os << "   Terminal Net Weighting:    " << tfStr[params.WSadjustment]
       << endl;
    if (params.WSadjustment) {
      os << "     Terminal Net Weighting Exponent:    "
         << params.termWeightModifier << endl;
    }
  }
  os << " Center-of-Gravity default cell location :  "
     << tfStr[params.cogLocation] << endl;
  os << " Capo Net Effect :  " << tfStr[params.capoNE] << endl;
  os << " Use Center of Gravity (COG) in SOR :  " << tfStr[!params.noCOG]
     << endl;
  os << " Use quadratic ini solution :  " << tfStr[params.useQuad] << endl;
  os << " Use linearization in QP :  " << tfStr[params.useLinearQP] << endl;
  os << " In SOR SkipNetsLargerThan : " << params.quadSkipNetsLargerThan << endl
     << endl;

  os << " Floorplanning Representation : " << params.FPrep << endl;
  os << " Treat macros as soft nodes : " << tfStr[params.allMacroSoft] << endl;
  if (params.allMacroSoft) {
    os << "  Maximum Aspect Ratio : " << params.maxAR << endl;
    os << "  Minimum Aspect Ratio : " << params.minAR << endl;
  }
#ifdef USEPATOMA
  os << " Use PATOMA : " << tfStr[params.usePatoma] << endl;
#endif
  os << " Anneal to cut lines : " << tfStr[params.annealToCutLines] << endl
     << endl;

  os << " Vertical AR factor for small partitioning: " << params.vertFactor
     << endl;
  os << " AR Stretch factor for large partitioning: " << params.ARStretchFactor
     << endl;

  os << endl << "----Bin Splitting Parameters----" << endl;
  os << " LookAhead       :  " << params.lookAhead << endl;
  os << params.splitterParams << endl;

  os << endl << "----Overlapping Parameters----" << endl;
  os << " Do overlapping:      " << tfStr[params.doOverlapping] << endl;
  if (params.doOverlapping) {
    if (params.startOverlappingLayer != 0)
      os << " StartOverlappingLayer:    " << params.startOverlappingLayer
         << endl;
    if (params.endOverlappingLayer != UINT_MAX)
      os << " EndOverlappingLayer:      " << params.endOverlappingLayer << endl;
  }
  os << " Honor Group Constraints:     " << tfStr[params.useGrpConstr] << endl;
  os << " Use Predetermined Cutlines:  "
     << tfStr[params.usePredeterminedCutlines] << endl;
  os << " Align Cutlines To Obstacles: " << tfStr[params.alignCutline2Obs]
     << endl;
  os << " Keep HGraph in memory during partitioning if its size is less than "
     << 100. * params.keepHGraphFraction << "% of the original" << endl;

  os << endl << "----Small Placer Parameters----" << endl;
  os << params.smplParams << endl;

  os << endl << "----Data Tracking Parameters----" << endl << endl;
  os << " Plot bins:           " << tfStr[params.plotBins] << endl;
  os << " Save layer BBoxes:   " << tfStr[params.saveLayerBBoxes] << endl;
  os << " Save bins:           " << tfStr[params.saveBins] << endl;
  os << " Save small pl. prob  " << tfStr[params.saveSmallPlProb] << endl;
  os << " Save bins in floorplan format:         "
     << tfStr[params.saveBinsFloorplan] << endl;
  os << " Save intermediate placements before partitioning  "
     << tfStr[params.prePartSavePl] << endl;

  os << endl << params.mlParams;

  return os;
}

CapoParameters::~CapoParameters() {}

void CapoParameters::printHelp() {
  cout << " ******* CAPO PARAMETERS HELP *********" << endl
       << "== Top-down parameters ==\n"
       << "  -verb  Verbosity\n"
       << "  -stopAtBins #\n"
       << "  -saveAtBins #\n"
       << "  -smallPartThreshold  <unsigned>\n"
       << "  -smallPlaceThreshold <unsigned>\n"
       << "  -smallSplitThreshold <unsigned>\n"
       << "  -replaceSmallBins <Never|AtTheEnd|AtEveryLayer>\n"
       << "  -useActualPinLocs\n"
       << "  -tdCongestionCtl\n"
       << "  -noWtCut\n"
       << "  -wtCutThresh <unsigned>\n"
       << "  -cogLocation\n"
       << "  -capoNE\n"
       << "  -noCOG\n"
       << "  -useQuad\n"
       << "  -useLinearQP\n"
       << "  -useHPWL\n"
       << "  -useFastSt\n"
#ifdef USEFLUTE
       << "  -useFLUTE\n"
#endif
       << "  -useMST\n"
       << "  -doFastPlaceMoves\n"
       << "  -prePartFix (fix cells before partitioning using intermediate "
          "placements)\n"
       << "  -prePartFixMargin <double %>\n"
       << "  -prePartFixLimit <double %>\n"
       << "  -prePartSavePl\n"
       << "  -useACG\n"
       << "  -tightenACG\n"
#ifdef USEHMETIS
       << "  -noHMetis\n"
#endif
       << "  -noThunder (no global placement refinement)\n"
       << "  -thunderMoveMultiplier <unsigned>\n"
       << "  -scampiThresh (scampi activation threshold)\n"
       << "  -useAnalytFP (use SOR for FP-ing selected instances)\n"
       << "  -useBloating\n"
       << "  -noWSadjustment\n"
       << "  -termWeightModifier <double> [2.0]\n"
       << "  -setLayerSeeds (e.g. 0:1:1)\n"
       << "  -quadSkipNetsLargerThan <unsigned>\n"
       << "  -selectiveFloorplanning <bool>\n"
       << "  -noClusterMacros <bool>\n"
       << "  -FPEvadeObstacles <bool>\n"
       << "  -FPEvadeObstaclesThresh <double 0-100% of bin Area> [0]\n"
       << "  -lookAheadFP <bool> \n"
       << "  -noMS   (no Mixed Size)\n\n"

       << "== Bin Splitting Parameters ==" << endl
       << "  -repart          <bool>\n"
       << "  -defaultTol      <double %>\n"
       << "  -minTol          <double %>\n"
       << "  -defaultRepartTol <double %>\n"
       << "  -minRepartTol    <double %>\n"
       << "  -tolCap          <double %>\n"
       << "  -constTol        <double %> [10]\n"
       << "  -uniformWS       <bool>\n"
       << "  -safeWS          <double %> [5]\n"
       << "  -minLocalWS      <double %> [2.5]\n"
       << "  -smallWS         <double %> [2]\n"
       << "  -congestionShifting <bool>\n"
       << "  -branchingFactor <unsigned> [2]\n"
       << "  -useNameHier     <delimiting char> [none]\n"
       << "  -lookAhead       <unsigned> [0]\n"
       << "  -noRowSplits     <bool>\n"
       << "  -noRepartSmallWS <bool>\n"
       << "  -noRepartShiftCL <bool>\n"
       << "  -useQuadCluster  <bool>\n"
       << "  -useCutOracle    <bool>\n"
       << "  -useOracleDirection <bool>\n"
       << "  -boost           <bool> [false]\n"
       << "   -boostLayer     <unsigned> [6]\n"
       << "   -boostFactor    <unsigned> [2]\n"
       << "  -noFeedback      <bool> [false]\n"
       << "  -feedback        <bool> [false]\n"
       << "   -termIter       <int> [3]\n"
       << "   -contMode       <int> [1]\n"
       << "   -acc            <bool>\n"
       << "  -useSeqPair      <bool> [true]  (use Sequence Pair representaion "
          "for floorplanning)\n"
       << "  -useBstarTree    <bool> [false] (use B*-Tree representation for "
          "floorplanning)\n"
       << "  -allMacroSoft    <bool> [false]\n"
       << "   -minAR          <double> [0.5] (sets minimum allowed aspect "
          "ratio for -allMacroSoft)\n"
       << "   -maxAR          <double> [2.0] (sets maximum allowed aspect "
          "ratio for -allMacroSoft)\n"
#ifdef USEPATOMA
       << "  -usePatoma       <bool> [false] (use PATOMA using of standard "
          "floorplanning)\n"
#endif
       << "  -selectSoft      <string> (the .blocks file determine which "
          "blocks are soft)\n"
       << "  -vertFactor      <double> [7.0]\n"
       << "     (aspect ratio threshold for switching a vertical cut to "
          "horizontal before small partitioning)\n"
       << "  -ARStretchFactor <double> [depends upon row spacing]\n"
       << "     (aspect ratio threshold for switching a vertical cut to "
          "horizontal before large partitioning)\n"
       << "  -annealToCutLines <bool> [false] tell the floorplanner to pack "
          "towards cut lines\n"
       << "  -termMergeThreshold <double %> [0]\n"
       << "     (if the percentage of terminals at the top level is greater "
          "than this, make a new subgraph before partitioning)\n\n"

       << "== Optional info-tracking parameters ==\n"
       << "  -plotBins       <bool> [false]\n"
       << "  -saveLayerBBoxes  <bool> [false]\n"
       << "  -saveBins       <bool> [false]\n"
       << "  -saveSmallPlProb  <bool> [false]\n"
       << "  -saveBinsFloorplan <bool> [false]\n\n"

       << "== Overlapping parameters ==\n"
       << "  -doOverlapping\n"
       << "  -startOverlappingLayer <unsigned>\n"
       << "  -endOverlappingLayer <unsigned>\n\n"

       << "== Constraints parameters ==\n"
       << "  -noRgnConstr\n"
       << "  -noGrpConstr\n"
       << "  -noPredeterminedCutlines\n"
       << "  -alignCutline2Obs\n"
       << "  -fillerCellsPartCaps\n\n"
       //        << "== Small Placer Parameters ==\n"
       //	  << smplParams

       << "== Meta-Option parameters ==\n"
       << "  -tryHarder  <unsigned>\n"
       << "  -competition  <bool>  [false]\n"
       << "  -faster    <unsigned>\n"

#ifdef USETRAFFIC
       << "== TRAFFIC parameters==\n"
       << "  -TRAFFIC\n     (enables use of TRAFFIC, off by default)\n"
       << "  -TRAFFICodds      <double> [50]\n     (chance of rotating any "
          "given block at random, 0-100)\n"
       << "  -TRAFFICthresh    <double> [10]\n     (threshold for grouping, "
          "0-100)\n"
#endif

       << endl;
}

// ---------------------------------------------------
// implementation of CapoSplitterParams

CapoSplitterParams::CapoSplitterParams(int argc, const char *argv[])
    : doRepartitioning(false),
      repartSmallWS(true),
      repartShiftCL(true),
      useQuadCluster(false),
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
      verb(argc, argv),
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
  BoolParam repart_("repart", argc, argv);
  BoolParam noRepartSmallWS_("noRepartSmallWS", argc, argv);
  BoolParam noRepartShiftCL_("noRepartShiftCL", argc, argv);
  BoolParam useQuadCluster_("useQuadCluster", argc, argv);
  BoolParam noRgnConstr_("noRgnConstr", argc, argv);
  BoolParam fillerCellsPartCaps_("fillerCellsPartCaps", argc, argv);
  StringParam useNameHier_("useNameHier", argc, argv);

  UnsignedParam tryharder_("tryHarder", argc, argv);
  UnsignedParam faster_("faster", argc, argv);
  BoolParam routable_("routable", argc, argv);
  BoolParam ROOSTER("ROOSTER", argc, argv);
  BoolParam FLROOSTER("FLROOSTER", argc, argv);
  DoubleParam ispd06("ispd06", argc, argv);
  DoubleParam densityTarget("densityTarget", argc, argv);

  // process meta-options first so they can be tweaked with other options
  if (ispd06.found()) {
    cout << "Meta-Option \"ispd06\" detected:" << endl;
    double target = ispd06;
    if (target > 1.) target *= 0.01;

    cout << "  - target density is " << 100. * target << "%" << endl << endl;

    minLocalWS = 1. - target;
    safeWS = 2. * minLocalWS;
  }

  if (densityTarget.found()) {
    double target = densityTarget;
    if (target > 1.) target *= 0.01;

    minLocalWS = 1. - target;
    safeWS = 2. * minLocalWS;
  }

  if (routable_.found() || ROOSTER.found() || FLROOSTER.found()) {
    // uniform tolerance on
    safeWS = -1;
    minLocalWS = -1;
    // congestion shifting on
    congestionShifting = true;
    if (routable_.found()) {
      cout << "Meta-Option \"routable\" detected:" << endl;
    }
    if (ROOSTER.found()) {
      cout << "Meta-Option \"ROOSTER\" detected:" << endl;
    }
    if (FLROOSTER.found()) {
      cout << "Meta-Option \"FLROOSTER\" detected:" << endl;
    }
    cout << "  - uniform partitioning tolerances enabled" << endl;
    cout << "  - congestion based WS allocation enabled" << endl;
  }

  if (tryharder_.found()) {
    if (tryharder_ >= 2) {
      numMLSets = max(tryharder_ + 1, 10u);
      cout << "Meta-Option \"tryHarder " << tryharder_
           << "\" detected:" << endl;
      cout << "  - branching factor set to " << numMLSets << endl;
    }
  }

  if (faster_.found()) {
    if (faster_ >= 2) {
      numMLSets = 1;
      cout << "Meta-Option \"faster " << faster_ << "\" detected:" << endl;
      cout << "  - branching factor set to " << numMLSets << endl;
    }
  }
  // end meta-options

  if (repart_.found()) doRepartitioning = true;

  if (noRepartSmallWS_.found()) repartSmallWS = false;

  if (noRepartShiftCL_.found()) repartShiftCL = false;

  if (useQuadCluster_.found()) useQuadCluster = true;

  if (noRgnConstr_.found()) useRgnConstr = false;

  if (fillerCellsPartCaps_.found()) fillerCellsPartCaps = true;

  if (useNameHier_.found()) {
    useNameHierarchy = true;
    delimiters = useNameHier_;
  }

  DoubleParam constTol_("constTol", argc, argv);
  DoubleParam defaultTol_("defaultTol", argc, argv);
  DoubleParam minTol_("minTol", argc, argv);
  DoubleParam defaultRepartTol_("defaultRepartTol", argc, argv);
  DoubleParam minRepartTol_("minRepartTol", argc, argv);
  DoubleParam tolCap_("tolCap", argc, argv);

  if (constTol_.found()) {
    useWSTolMethod = false;
    defaultTolerance = constTol_ / 100.;
  }

  if (defaultTol_.found()) defaultTolerance = defaultTol_ / 100.;

  abkfatal(defaultTolerance >= 0.,
           "The default tolerance must be non-negative");

  if (minTol_.found()) minTolerance = minTol_ / 100.;

  abkfatal(minTolerance >= 0., "The minimum tolerance must be non-negative");
  abkfatal(defaultTolerance >= minTolerance,
           "Default tolerance must be at least as much as minimum tolerance");

  if (defaultRepartTol_.found())
    defaultRepartTolerance = defaultRepartTol_ / 100.;

  abkfatal(defaultRepartTolerance >= 0.,
           "The default repartitioning tolerance must be non-negative");

  if (minRepartTol_.found()) minRepartTolerance = minRepartTol_ / 100.;

  abkfatal(minRepartTolerance >= 0.,
           "The minimum repartitioning tolerance must be non-negative");
  abkfatal(defaultRepartTolerance >= minRepartTolerance,
           "Default repartitining tolerance must be at least as much as "
           "minimum repartitioning tolerance");

  if (tolCap_.found()) tolCap = tolCap_ / 100.;

  abkfatal(tolCap >= 0., "Tolerance cap must be non-negative");

  DoubleParam safeWS_("safeWS", argc, argv);
  DoubleParam minLocalWS_("minLocalWS", argc, argv);
  BoolParam uniformWS_("uniformWS", argc, argv);
  DoubleParam smallWS_("smallWS", argc, argv);
  BoolParam congestionShifting_("congestionShifting", argc, argv);
  BoolParam alternateCongShifting("alternateCongShifting", argc, argv);
  DoubleParam congExponent_("congExponent", argc, argv);

  if (safeWS_.found()) {
    safeWS = safeWS_ / 100.;
    abkwarn(safeWS >= 0.,
            "Negative safe white space disables non uniform white space");
  }

  if (minLocalWS_.found()) {
    minLocalWS = minLocalWS_ / 100.;
    abkwarn(
        minLocalWS >= 0.,
        "Negative min local white space disables minimum whitespace checking");
  }

  abkfatal(safeWS < 0 || minLocalWS <= safeWS,
           "Safe white space must be larger than min local white space");

  if (uniformWS_.found()) {
    abkwarn((safeWS < 0 || !safeWS_.found()) &&
                (minLocalWS < 0 || !minLocalWS_.found()),
            "Conflicting options found. uniformWS will override safeWS and "
            "minLocalWS");
    safeWS = -1;
    minLocalWS = -1;
  }

  if (congestionShifting_.found()) {
    congestionShifting = true;
  }

  if (alternateCongShifting.found()) {
    alwaysCongShift = false;
  }

  if (congExponent_.found()) {
    congExponent = congExponent_;
  }

  if (smallWS_.found()) {
    smallWS = smallWS_ / 100.;
    abkfatal(smallWS >= 0., "Small white space must be non-negative");
  }

  if (safeWS >= 0) smallWS = min(smallWS, safeWS);
  if (minLocalWS >= 0) smallWS = min(smallWS, minLocalWS);

  UnsignedParam bFactor_("branchingFactor", argc, argv);
  DoubleParam bFactorAdjustCutoff_("branchingFactorAdjustCutoff", argc, argv);

  if (bFactor_.found()) numMLSets = bFactor_;

  if (bFactorAdjustCutoff_.found()) bFactorAdjustCutoff = bFactorAdjustCutoff_;

  DoubleParam extraBranchCutoff_("extraBranchCutoff", argc, argv);
  DoubleParam termMergeThreshold_("termMergeThreshold", argc, argv);
  BoolParam newFuzzy_("newFuzzy", argc, argv);
  BoolParam newFuzzyOnlyRepart_("newFuzzyOnlyRepart", argc, argv);

  if (extraBranchCutoff_.found()) extraBranchCutoff = extraBranchCutoff_;

  if (termMergeThreshold_.found())
    termMergeThreshold = termMergeThreshold_ / 100.;

  if (newFuzzy_.found()) newFuzzy = newFuzzy_;

  if (newFuzzyOnlyRepart_.found()) {
    newFuzzy = true;
    newFuzzyOnlyRepart = newFuzzyOnlyRepart_;
  }

  BoolParam useCutOracle_("useCutOracle", argc, argv);
  BoolParam useOracleDirection_("useOracleDirection", argc, argv);
  if (useCutOracle_.found()) useCutOracle = true;

  if (useOracleDirection_.found()) useOracleDirection = true;

  BoolParam useACG_("useACG", argc, argv);
  BoolParam tightenACG_("tightenACG", argc, argv);
  if (useACG_.found()) {
    useACG = true;
  }
  if (tightenACG_.found()) tightenACG = true;

  BoolParam WSOracle_("WSOracle", argc, argv);
  if (WSOracle_.found()) WSOracle = true;

  DoubleParam WSOracleTol_("WSOracleTol", argc, argv);
  if (WSOracleTol_.found()) {
    WSOracle = true;
    WSOracleTol = WSOracleTol_ / 100.;
  }

  BoolParam eco_("ECO", argc, argv);
  if (eco_.found()) {
    eco = true;
  }

  BoolParam ecoShift_("ecoShift", argc, argv);
  if (ecoShift_.found()) {
    eco = true;
    ecoShift = true;
  }

#ifdef USEFLOW
  BoolParam flowECO_("flowECO", argc, argv);
  if (flowECO_.found()) {
    eco = true;
    flowECO = true;
  }
#endif

  if (eco && ispd06.found()) {
    ecoShift = true;
  }

  DoubleParam ecoOverfillPct_("ecoOverfillPct", argc, argv);
  if (ecoOverfillPct_.found()) {
    ecoOverfillPct = 0.01 * ecoOverfillPct_;
  }

  DoubleParam ecoImprovePct_("ecoImprovePct", argc, argv);
  if (ecoImprovePct_.found()) {
    ecoImprovePct = 0.01 * ecoImprovePct_;
  }

  abkfatal(!(WSOracle && eco),
           "WSOracle and ECO mode are currently incompatible together");

#ifdef USEHMETIS
  BoolParam noHMetis_("noHMetis", argc, argv);
  if (noHMetis_.found()) useHMetis = false;

  DoubleParam hmetisProblemThresh_("hmetisProblemThresh", argc, argv);
  if (hmetisProblemThresh_.found())
    hmetisProblemThresh = hmetisProblemThresh_ / 100.;
#endif

  BoolParam useFMPartPlus_("useFMPartPlus", argc, argv);
  if (useFMPartPlus_.found()) useFMPartPlus = true;

  // <aaronnn> adjust WS after partitioning by utilization density
  BoolParam adjustWSByDensity_("adjustWSByDensity", argc, argv);
  if (adjustWSByDensity_.found()) adjustWSByDensity = true;

  UnsignedParam cellGroups("numCellGroups", argc, argv);
  if (cellGroups.found()) {
    numCellGroups = cellGroups;
  }

  UnsignedParam cellGroupsCutoff_("cellGroupsCutoff", argc, argv);
  if (cellGroupsCutoff_.found()) {
    cellGroupsCutoff = cellGroupsCutoff_;
  }

  DoubleParam netDegreeWeightExponent_("netDegreeWeightExponent", argc, argv);
  if (netDegreeWeightExponent_.found()) {
    netDegreeWeighting = true;
    netDegreeWeightExponent = netDegreeWeightExponent_;
  }
}

CapoSplitterParams::CapoSplitterParams(const CapoSplitterParams &srcParams)
    : doRepartitioning(srcParams.doRepartitioning),
      repartSmallWS(srcParams.repartSmallWS),
      repartShiftCL(srcParams.repartShiftCL),
      useQuadCluster(srcParams.useQuadCluster),
      useWSTolMethod(srcParams.useWSTolMethod),
      defaultTolerance(srcParams.defaultTolerance),
      defaultRepartTolerance(srcParams.defaultRepartTolerance),
      minTolerance(srcParams.minTolerance),
      minRepartTolerance(srcParams.minRepartTolerance),
      tolCap(srcParams.tolCap),
      safeWS(srcParams.safeWS),
      minLocalWS(srcParams.minLocalWS),
      smallWS(srcParams.smallWS),
      congestionShifting(srcParams.congestionShifting),
      alwaysCongShift(srcParams.alwaysCongShift),
      congExponent(srcParams.congExponent),
      WSOracle(srcParams.WSOracle),
      WSOracleTol(srcParams.WSOracleTol),
      eco(srcParams.eco),
      ecoShift(srcParams.ecoShift),
#ifdef USEFLOW
      flowECO(srcParams.flowECO),
#endif
      ecoOverfillPct(srcParams.ecoOverfillPct),
      ecoImprovePct(srcParams.ecoImprovePct),
      numMLSets(srcParams.numMLSets),
      useNameHierarchy(srcParams.useNameHierarchy),
      delimiters(srcParams.delimiters),
      useRgnConstr(srcParams.useRgnConstr),
      adjustWSByDensity(srcParams.adjustWSByDensity),
      fillerCellsPartCaps(srcParams.fillerCellsPartCaps),
      verb(srcParams.verb),
      bFactorAdjustCutoff(srcParams.bFactorAdjustCutoff),
      extraBranchCutoff(srcParams.extraBranchCutoff) /*,
       savePartSolnWhenUndoing(srcParams.savePartSolnWhenUndoing)*/
      ,
      termMergeThreshold(srcParams.termMergeThreshold),
      newFuzzy(srcParams.newFuzzy),
      newFuzzyOnlyRepart(srcParams.newFuzzyOnlyRepart),
      useCutOracle(srcParams.useCutOracle),
      useOracleDirection(srcParams.useOracleDirection),
      useACG(srcParams.useACG),
      tightenACG(srcParams.tightenACG),
#ifdef USEHMETIS
      useHMetis(srcParams.useHMetis),
      hmetisProblemThresh(srcParams.hmetisProblemThresh),
#endif
      useFMPartPlus(srcParams.useFMPartPlus),
      numCellGroups(srcParams.numCellGroups),
      cellGroupsCutoff(srcParams.cellGroupsCutoff),
      netDegreeWeighting(srcParams.netDegreeWeighting),
      netDegreeWeightExponent(srcParams.netDegreeWeightExponent) {
}

ostream &operator<<(ostream &os, const CapoSplitterParams &params) {
  const char *tfStr[2];
  tfStr[0] = "false";
  tfStr[1] = "true";

  os << "---Capo Splitter Parameters---" << endl;
  os << "  doRepartitioning   " << tfStr[params.doRepartitioning] << endl;
  os << "  noRepartSmallWS    " << tfStr[!params.repartSmallWS] << endl;
  os << "  noRepartShiftCL    " << tfStr[!params.repartShiftCL] << endl;
  os << "  Whitespace allocation: " << endl;
  os << "    congestion whitespace: " << tfStr[params.congestionShifting]
     << endl;
  os << "    uniform whitespace: "
     << tfStr[params.safeWS < 0. && params.minLocalWS < 0.] << endl;
  if (params.safeWS >= 0.)
    os << "    safeWS     " << params.safeWS * 100. << "%" << endl;
  if (params.minLocalWS >= 0.)
    os << "    minLocalWS " << params.minLocalWS * 100. << "%" << endl;
  os << "  adjustWSByDensity    " << tfStr[params.adjustWSByDensity] << endl;

  if (params.useWSTolMethod)
    os << "  h-cuts use WS Tolerance Method     " << endl;
  else
    os << "  h-cuts use constant tolerance " << params.defaultTolerance * 100.
       << "%" << endl;

  if (params.useWSTolMethod)
    os << "  v-cuts use WS Tolerance Method     " << endl;
  else
    os << "  v-cuts use constant tolerance " << params.defaultTolerance * 100.
       << "%" << endl;

  os << "  Default Tolerance " << params.defaultTolerance * 100. << "%" << endl;
  os << "  Min Tolerance     " << params.minTolerance * 100. << "%" << endl;
  os << "  Default Repart Tolerance " << params.defaultRepartTolerance * 100.
     << "%" << endl;
  os << "  Min Repart Tolerance " << params.minRepartTolerance * 100. << "%"
     << endl;
  os << "  Tolerance Cap " << params.tolCap * 100. << "%" << endl;

  os << "  WSOracle " << tfStr[params.WSOracle] << endl;
  if (params.WSOracle)
    os << "    WSOracle Tolerance " << params.WSOracleTol * 100. << "%" << endl;

  os << "  ECO mode " << tfStr[params.eco] << endl;
  if (params.eco) {
    os << "    ECO WS " << 100. * params.ecoOverfillPct << "%" << endl;
    os << "    ECO Improvement Threshold " << 100. * params.ecoImprovePct << "%"
       << endl;
  }

  os << "  use name hierarchy    " << tfStr[params.useNameHierarchy] << endl;
  if (params.useNameHierarchy)
    os << "   delimiter          " << params.delimiters << endl;
  os << "  Respect Region Utilization Constraints : "
     << tfStr[params.useRgnConstr] << endl;
  os << "  With Filler Cells Put Std-Cells In Center : "
     << tfStr[params.fillerCellsPartCaps] << endl;
  os << "  Use Cut Position Oracle :  " << tfStr[params.useCutOracle] << endl;
  os << "  Use Cut Direction Oracle : " << tfStr[params.useOracleDirection]
     << endl;
  os << "  Use Analytical-Constraint-Generation : " << tfStr[params.useACG]
     << endl;
  os << "  Tighten Analytical-Constraint-Generation : "
     << tfStr[params.tightenACG] << endl;
#ifdef USEHMETIS
  os << "  Using HMetis : " << tfStr[params.useHMetis] << endl;
  if (params.useHMetis)
    os << "    HMetis Problem Detection Threshold : "
       << params.hmetisProblemThresh * 100 << "%" << endl;
#endif
  os << "  termMergeThreshold " << params.termMergeThreshold << endl;

  return os;
}
