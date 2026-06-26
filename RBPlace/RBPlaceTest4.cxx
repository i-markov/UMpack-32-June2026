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

// Row-based legalization and plotting driver.
// It exposes the legalization, snapping, tethering, plotting, and save
// options used by the row-based placement flow on a single test executable.

#include <iostream>

#include "ABKCommon/abkassert.h"
#include "ABKCommon/abktimer.h"
#include "ABKCommon/paramproc.h"
#include "RBPlacePlot.h"
#include "RBPlacement.h"

using std::cout;
using std::endl;
using std::vector;

int main(int argc, const char *argv[]) {
  NoParams noParams(argc, argv);
  BoolParam helpRequest("help", argc, argv);
  BoolParam helpRequest1("h", argc, argv);
  BoolParam remOverlaps("legal", argc, argv);
  BoolParam checkLegal("checkLegal", argc, argv);
  BoolParam markMacros("markMacros", argc, argv);
  BoolParam snapToRows("snapToRows", argc, argv);
  BoolParam snapToSites("snapToSites", argc, argv);
  BoolParam stdNodeHeights("stdNodeHeights", argc, argv);

  RBPlacement::Parameters rbParams(argc, argv);
  DoubleParam addDummy("addDummy", argc, argv);
  DoubleParam tetherCells("tetherCells", argc, argv);
  BoolParam takeTetherConstrFrmFile("takeTetherConstrFrmFile", argc, argv);
  BoolParam takeTetherNetConstrFrmFile("takeTetherNetConstrFrmFile", argc,
                                       argv);
  DoubleParam rgnSizePercent("rgnSizePercent", argc, argv);
  DoubleParam tetherNewAR("tetherNewAR", argc, argv);
  DoubleParam tetherNewWS("tetherNewWS", argc, argv);
  StringParam saveAsNodes("saveAsNodes", argc, argv);
  StringParam saveAsCplace("saveAsCplace", argc, argv);
  StringParam saveAsPlato("saveAsPlato", argc, argv);
  StringParam saveLEFDEF("saveLEFDEF", argc, argv);
  BoolParam markMacrosAsBlocks("markMacrosAsBlocks", argc, argv);
  BoolParam remSitesMacro("remSitesMacro", argc, argv);
  BoolParam fixMacros("fixMacros", argc, argv);
  DoubleParam remSitesCong("remSitesCong", argc, argv);
  DoubleParam addNodesCong("addNodesCong", argc, argv);
  StringParam saveAsNodesFloorplan("saveAsNodesFloorplan", argc, argv);
  StringParam saveMacrosAsNodesFloorplan("saveMacrosAsNodesFloorplan", argc,
                                         argv);
  UnsignedParam useSingleNet("useSingleNet", argc, argv);
  DoubleParam maxCellHeight("maxCellHeight", argc, argv);
  DoubleParam maxCellWidth("maxCellWidth", argc, argv);
  StringParam saveAsNodesShredHW("saveAsNodesShredHW", argc, argv);
  StringParam saveAsNodesWPinAssign("saveAsNodesWPinAssign", argc, argv);
  StringParam savePlacementUnShred("savePlacementUnShred", argc, argv);

  StringParam savePl("savePl", argc, argv);
  StringParam savePlNoDummy("savePlNoDummy", argc, argv);
  StringParam plotNets("plotNets", argc, argv);
  StringParam plotNodes("plotNodes", argc, argv);
  StringParam plotNodesWNoFiller("plotNodesWNoFiller", argc, argv);
  StringParam plotNodesAndNets("plot", argc, argv);
  StringParam plotNodesWNames("plotNodesWNames", argc, argv);
  StringParam plotNodesAndNetsWNames("plotWNames", argc, argv);
  StringParam plotSites("plotSites", argc, argv);
  StringParam plotNodesWSites("plotNodesWSites", argc, argv);
  StringParam plotRows("plotRows", argc, argv);
  StringParam plotDummy("plotDummy", argc, argv);
  StringParam plotWSMap("plotWSMap", argc, argv);
  StringParam plotFSFP("plotFSFP", argc, argv);  // Free-shape Floorplan

  BoolParam spaceCellsEqually("spaceCellsEqually", argc, argv);
  BoolParam spaceCellsWCongInfo("spaceCellsWCongInfo", argc, argv);
  BoolParam shiftCellsLeft("shiftCellsLeft", argc, argv);

  BoolParam assignPinsToCells("assignPinsToCells", argc, argv);

  if (noParams.found() || helpRequest.found() || helpRequest1.found()) {
    cout << "  Use '-help' or '-f filename.aux' " << endl;
    cout << "  Other options: \n"
         << "   -remSitesMacro   (remove Sites below Macros)\n"
         << "   -remSitesCong double (remove % of whitespace from congested "
            "regions)\n"
         << "   -addNodesCong double (attach % of whitespace dummy nodes to "
            "congested nodes)\n"
         << "   -stdNodeHeights  (set all node heights to row height)\n"
         << "   -fixMacros   (mark Macros as terminals while saving)\n"
         << "   -tetherCells double (add tethering nets to fract % of cells to "
            "achieve stability)\n"
         << "   -rgnSizePercent double (constraining region size for each cell "
            "during tethering. % of layout)\n"
         << "   -takeTetherConstrFrmFile (Take part of the nodes to tether "
            "from constraints.nodes input file\n"
         << "   -takeTetherNetConstrFrmFile (Don't tether nodes connected to "
            "nets specified in constraints.nets input file\n"
         << "   -tetherNewAR  (Tether for a new reshaped layout)\n"
         << "   -tetherNewWS  (Tether for a new rescaled layout)\n"
         << "   -addDummy double (add dummy filler cells. keep min % of "
            "whitespace)\n"
         << "   -saveAsNodes     base_filename\n"
         << "   -saveLEFDEF      base_filename\n"
         << "   -saveAsCplace    base_filename\n"
         << "   -saveAsPlato     base_filename\n"
         << "   -markMacrosAsBlocks (for saveLEFDEF)\n"
         << "   -assignPinsToCells \n"
         << "   -savePl  filename\n"
         << "   -savePlNoDummy  filename  (savePl with no dummy objects)\n"
         << "   -saveAsNodesWPinAssign     base_filename\n"
         << "   -saveAsNodesFloorplan base_filename\n"
         << "   -saveMacrosAsNodesFloorplan base_filename\n"
         << "   -saveAsNodesShredHW  base_filename\n"
         << "   -useSingleNet   unsigned (use weighted single net construction "
            "while shredding)\n"
         << "   -maxCellWidth   double (width of shredded cells)\n"
         << "   -maxCellHeight  double (height of shredded cells)\n"
         << "   -savePlacementUnShred  filename\n"
         << "   -xmin,-xmax,-ymin,-ymax  (for plotting range. Default means "
            "all)\n"
         << "   -plotNets        base_fileName\n"
         << "   -plotNodes       base_fileName\n"
         << "   -plotNodesWNoFiller base_fileName \n"
         << "   -plotNodesWNames base_fileName\n"
         << "   -plot            base_fileName   (plots nodes and nets)\n"
         << "   -plotWNames      base_fileName   (the above + node names)\n"
         << "   -plotSites	  base_fileName   (plots the site map)\n"
         << "   -plotNodesWSites base_fileName   (the above + nodes)\n"
         << "   -plotRows        base_fileName   (plot the rows)\n"
         << "   -plotDummy       base_fileName   (plot only dummy nodes+nets)\n"
         << "   -plotWSMap       base_fileName   (plot the whitespace "
            "histogram)\n"
         << "   -plotFSFP        base_fileName   (plot the free-shape "
            "floorplan)\n"
         << "   -checkLegal (Use to check for out of core cells / overlaps)\n"
         << "   -printOverlaps (Print out the bounding boxes of all the "
            "overlaps)\n"
         << "   -legal  (Use to remove any overlaps in existing placement)\n"
         << "   -snapToRows  (Use to snap cells to closest rows. e.g. for "
            "Kraftwerk Placements )\n"
         << "   -snapToSites (Use to snap cells to sites in existing "
            "placement)\n"
         << "   -snapMacros  (Snap macros to rows and sites during macro "
            "overlap removal)\n"
         << "   -spaceCellsEqually   (for each subrow, space cells equally)\n"
         << "   -spaceCellsWCongInfo (for each subrow, space cells w "
            "congestion info)\n"
         << "   -shiftCellsLeft      (pack Cells to the left of each subrow)\n"
         << endl;
    return 0;
  }
  cout << "RBPlacement Parameters" << endl;
  cout << rbParams << endl;

  StringParam auxFile("f", argc, argv);
  abkfatal(auxFile.found(), "must have -f <auxfilename>");
  RBPlacement rbplace(auxFile, rbParams);

  if (assignPinsToCells.found()) {
    rbplace.assignPinsToCells(true);
  }

  if (snapToRows.found()) {
    rbplace.alignCellsToRows();
    rbplace.snapCellsInSites();
  }

  if (stdNodeHeights.found()) {
    double rowHeight = (rbplace.coreRowsBegin())->getHeight();
    const_cast<HGraphWDimensions &>(rbplace.getHGraph()).setNodeDims(rowHeight);
  }

  if (spaceCellsWCongInfo.found()) rbplace.spaceCellsWCongInfoInRows1();

  if (spaceCellsEqually.found()) rbplace.spaceCellsEquallyInRows();

  if (shiftCellsLeft.found()) rbplace.shiftCellsLeft();

  if (remSitesCong.found()) rbplace.removeSitesFromCongestedRgn(remSitesCong);

  MaxMem maxMem;

#ifdef USEFLOW
  BoolParam flow("flow", argc, argv);
  BoolParam xFlow("xFlow", argc, argv);
  BoolParam yFlow("yFlow", argc, argv);

  bool xflows = false;
  bool yflows = false;

  if (flow.found()) {
    xflows = true;
    yflows = true;
  }

  if (xFlow.found()) {
    xflows = true;
  }

  if (yFlow.found()) {
    yflows = true;
  }

  if (xflows || yflows) {
    std::pair<double, double> HPWL = rbplace.evalXYHPWL(false);
    cout << " Center-to-center HalfPerim WL: " << HPWL.first + HPWL.second
         << " (" << HPWL.first << " x , " << HPWL.second << " y )" << endl;
    HPWL = rbplace.evalXYHPWL(true);
    cout << " Pin-to-pin HalfPerim WL      : " << HPWL.first + HPWL.second
         << " (" << HPWL.first << " x , " << HPWL.second << " y )" << endl;

    Timer flowTime;
    if (rbplace.getNumMacros() > 0) {
      rbplace.removeSitesBelowMacros();
      rbplace.alignCellsToRows();
    }

    rbplace.doFlow(NULL, xflows, yflows, &maxMem);

    if (rbplace.getNumMacros() > 0) {
      rbplace.reinstateCoreRowsAndRepopulate();
    }
    flowTime.stop();
    cout << "Flows took " << flowTime.getUserTime() << endl;
  }
#endif

#ifdef USEFLOW
  BoolParam unconstrX("unconstrainedX", argc, argv);
  BoolParam unconstrY("unconstrainedY", argc, argv);

  if (unconstrX.found() || unconstrY.found()) {
    vector<unsigned> cellIds;
    vector<BBox> bboxes;

    const BBox &theBBox = rbplace.getBBox(false);

    PlacementWOrient temp = rbplace.getPlacement();

    for (unsigned i = rbplace.getHGraph().getNumTerminals();
         i < rbplace.getHGraph().getNumNodes(); ++i) {
      if (!rbplace.isFixed(i)) {
        cellIds.push_back(i);
        bboxes.push_back(theBBox);
      }
    }

    if (unconstrX.found()) {
      rbplace.doUnconstrainedXFlow(rbplace.getPlacement(), rbplace.getHGraph(),
                                   cellIds, bboxes, temp, &maxMem);
    }

    if (unconstrY.found()) {
      rbplace.doUnconstrainedYFlow(rbplace.getPlacement(), rbplace.getHGraph(),
                                   cellIds, bboxes, temp, &maxMem);
    }

    rbplace.updatePlacementWOri(temp);
  }
#endif

  if (remOverlaps.found()) {
    if (rbplace.getNumMacros() > 0) {
      rbplace.removeSitesBelowMacros();
      rbplace.alignCellsToRows();

      BoolParam snap("snapMacros", argc, argv);
      bool snapMacrosX = snap.found();
      bool snapMacrosY = snap.found();
      rbplace.removeMacroOverlaps(&maxMem, snapMacrosX, snapMacrosY);
    }
    rbplace.snapCellsInSites();
    rbplace.royjRemOverlaps(&maxMem);
    if (rbplace.getNumMacros() > 0) {
      rbplace.reinstateCoreRowsAndRepopulate();
    }
  }

  if (snapToSites.found()) rbplace.snapCellsInSites();

  if (checkLegal.found()) {
    unsigned numCellsNotPlaced = rbplace.numCellsNotPlaced();
    unsigned numOverlaps = rbplace.getNumOverlaps();

    rbplace.checkOutOfCoreCells();
    cout << " RBPlacement: Overlap Statistics\n";
    BoolParam printOverlaps("printOverlaps", argc, argv);
    rbplace.calcOverlapGeneric(OnlyMacro);
    double overlap = rbplace.calcOverlapGeneric(Movable);
    rbplace.calcOverlapGeneric(OnlyFixed);
    rbplace.calcOverlapGeneric(All, printOverlaps);

    if (numCellsNotPlaced > 0 || numOverlaps > 0 ||
        greaterThanDouble(overlap, 0.))
      cout << " Output of legality checker : === FALSE ===" << endl;
    else
      cout << " Output of legality checker : === TRUE ===" << endl;
    cout << endl;
  }

  if (remSitesMacro) rbplace.removeSitesBelowMacros();

  if (saveLEFDEF.found()) rbplace.saveLEFDEF(saveLEFDEF, markMacrosAsBlocks);

  if (saveAsNodes.found()) {
    if (addDummy.found()) {
      rbplace.saveAsNodesNetsWtsWDummy(saveAsNodes, addDummy);
    } else if (addNodesCong.found()) {
      rbplace.saveAsNodesNetsWtsWCongInfo(saveAsNodes, addNodesCong);
    } else if (tetherCells.found()) {
      double regionSizePercent = 0.005;
      double newAR = -1;
      double newWS = -1;
      if (rgnSizePercent.found()) regionSizePercent = rgnSizePercent;
      if (tetherNewAR.found()) newAR = tetherNewAR;
      if (tetherNewWS.found()) newWS = tetherNewWS;

      rbplace.saveAsNodesNetsWtsWTether(
          saveAsNodes, tetherCells, regionSizePercent, takeTetherConstrFrmFile,
          takeTetherNetConstrFrmFile, newAR, newWS);
    } else {
      if (fixMacros.found())
        rbplace.saveAsNodesNetsWts(saveAsNodes, fixMacros);
      else
        rbplace.saveAsNodesNetsWts(saveAsNodes);
    }
  }

  if (saveAsNodesWPinAssign.found()) {
    rbplace.assignPinsToCells();
    rbplace.saveAsNodesNetsWts(saveAsNodesWPinAssign);
  }

  if (saveAsNodesFloorplan.found())
    rbplace.saveAsNodesNetsPlFloorplan(saveAsNodesFloorplan);

  if (saveMacrosAsNodesFloorplan.found())
    rbplace.saveMacrosAsNodesNetsPlFloorplan(saveMacrosAsNodesFloorplan);

  if (saveAsNodesShredHW.found()) {
    unsigned singleNetWt = 0;
    double cellWidth = 0;
    double cellHeight = 0;
    if (useSingleNet.found()) singleNetWt = useSingleNet;
    if (maxCellWidth.found()) cellWidth = maxCellWidth;
    if (maxCellHeight.found()) cellHeight = maxCellHeight;

    rbplace.saveAsNodesNetsWtsShredHW(saveAsNodesShredHW, cellHeight, cellWidth,
                                      singleNetWt);
  }

  if (savePlacementUnShred.found())
    rbplace.savePlacementUnShredHW(savePlacementUnShred);

  if (savePl.found()) rbplace.savePlacement(savePl);

  if (savePlNoDummy.found()) rbplace.savePlNoDummy(savePlNoDummy);

  if (saveAsCplace.found()) rbplace.saveAsCplace(saveAsCplace);

  if (saveAsPlato.found()) rbplace.saveAsPlato(saveAsPlato);

  Plotters::RBPlacePlotter::Parameters allPlotParams(argc, argv, false);
  if (plotNets.found()) {
    Plotters::RBPlacePlotter::Parameters plotParams(allPlotParams);
    plotParams.plotNets = true;
    rbplace.saveAsPlot(plotParams, plotNets);
  }

  if (plotNodes.found()) {
    Plotters::RBPlacePlotter::Parameters plotParams(allPlotParams);
    plotParams.plotNodes = true;
    rbplace.saveAsPlot(plotParams, plotNodes);
  }

  if (plotNodesWNoFiller.found()) {
    Plotters::RBPlacePlotter::Parameters plotParams(allPlotParams);
    plotParams.plotNodes = true;
    plotParams.plotFillerCells = false;
    rbplace.saveAsPlot(plotParams, plotNodesWNoFiller);
  }

  if (plotNodesAndNets.found()) {
    Plotters::RBPlacePlotter::Parameters plotParams(allPlotParams);
    plotParams.plotNodes = true;
    plotParams.plotNets = true;
    rbplace.saveAsPlot(plotParams, plotNodesAndNets);
  }

  if (plotNodesWNames.found()) {
    Plotters::RBPlacePlotter::Parameters plotParams(allPlotParams);
    plotParams.plotNodes = true;
    plotParams.plotNodesNames = true;
    rbplace.saveAsPlot(plotParams, plotNodesWNames);
  }

  if (plotNodesAndNetsWNames.found()) {
    Plotters::RBPlacePlotter::Parameters plotParams(allPlotParams);
    plotParams.plotNodes = true;
    plotParams.plotNodesNames = true;
    plotParams.plotNets = true;
    rbplace.saveAsPlot(plotParams, plotNodesAndNetsWNames);
  }

  if (plotSites.found()) {
    Plotters::RBPlacePlotter::Parameters plotParams(allPlotParams);
    plotParams.plotSiteMap = true;
    rbplace.saveAsPlot(plotParams, plotSites);
  }

  if (plotNodesWSites.found()) {
    Plotters::RBPlacePlotter::Parameters plotParams(allPlotParams);
    plotParams.plotSiteMap = true;
    plotParams.plotNodes = true;
    rbplace.saveAsPlot(plotParams, plotNodesWSites);
  }

  if (plotRows.found()) {
    Plotters::RBPlacePlotter::Parameters plotParams(allPlotParams);
    plotParams.plotRows = true;
    rbplace.saveAsPlot(plotParams, plotRows);
  }

  if (plotDummy.found()) {
    Plotters::RBPlacePlotter::Parameters plotParams(allPlotParams);
    plotParams.plotNodes = true;
    plotParams.plotNets = true;
    rbplace.plotDummy(plotParams, plotDummy);
  }

  if (plotWSMap.found()) rbplace.plotWSHist(plotWSMap);

  if (plotFSFP.found())
    Plotters::RBPlacePlotter::plotFreeshapeFloorplan(
        std::string(static_cast<const char *>(plotFSFP)), rbplace, argc, argv);

  std::pair<double, double> HPWL = rbplace.evalXYHPWL(false);
  cout << " Center-to-center HalfPerim WL: " << HPWL.first + HPWL.second << " ("
       << HPWL.first << " x , " << HPWL.second << " y )" << endl;
  HPWL = rbplace.evalXYHPWL(true);
  cout << " Pin-to-pin HalfPerim WL      : " << HPWL.first + HPWL.second << " ("
       << HPWL.first << " x , " << HPWL.second << " y )" << endl;

  // cout<<endl<<endl<<auxFile<<" numMacros is :
  // "<<rbplace.getNumMacros()<<endl<<endl;

  return 0;
}
