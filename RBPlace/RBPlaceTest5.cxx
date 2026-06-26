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
// It exercises the same placement utilities as Test4 but with a smaller
// option set for quick smoke testing of legalization and export paths.

#include <iostream>

#include "ABKCommon/abkassert.h"
#include "ABKCommon/paramproc.h"
#include "RBPlacePlot.h"
#include "RBPlacement.h"

using std::cout;
using std::endl;

int main(int argc, const char *argv[]) {
  NoParams noParams(argc, argv);
  BoolParam helpRequest("help", argc, argv);
  BoolParam helpRequest1("h", argc, argv);
  BoolParam remOverlaps("legal", argc, argv);
  BoolParam snapToSites("snapToSites", argc, argv);
  BoolParam stdNodeHeights("stdNodeHeights", argc, argv);

  RBPlacement::Parameters rbParams(argc, argv);
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
  StringParam saveAsNodesShredHW("saveAsNodesShredHW", argc, argv);
  StringParam saveAsNodesWPinAssign("saveAsNodesWPinAssign", argc, argv);
  StringParam savePlacementUnShred("savePlacementUnShred", argc, argv);
  StringParam savePl("savePl", argc, argv);
  StringParam savePlNoDummy("savePlNoDummy", argc, argv);
  StringParam plotNets("plotNets", argc, argv);
  StringParam plotNodes("plotNodes", argc, argv);
  StringParam plotNodesAndNets("plot", argc, argv);
  StringParam plotNodesWNames("plotNodesWNames", argc, argv);
  StringParam plotNodesAndNetsWNames("plotWNames", argc, argv);
  StringParam plotSites("plotSites", argc, argv);
  StringParam plotNodesWSites("plotNodesWSites", argc, argv);
  StringParam plotRows("plotRows", argc, argv);

  BoolParam spaceCellsEqually("spaceCellsEqually", argc, argv);
  BoolParam spaceCellsWCongInfo("spaceCellsWCongInfo", argc, argv);
  BoolParam shiftCellsLeft("shiftCellsLeft", argc, argv);

  cout
      << "/********************************************************************"
         "******\n"
      << "***\n"
      << "*** Copyright (c) 1995-2000 Regents of the University of "
         "California,\n"
      << "***               Andrew E. Caldwell, Andrew B. Kahng and Igor L. "
         "Markov\n"
      << "*** Copyright (c) 2000-2006 Regents of the University of Michigan,\n"
      << "***               Saurabh N. Adya, Jarrod A. Roy and Igor L. Markov\n"
      << "***\n"
      << "***  Contact author(s): abk@cs.ucsd.edu, igor.markov1@gmail.com\n"
      << "***  Original Affiliation:   UCLA, Computer Science Department,\n"
      << "***                          Los Angeles, CA 90095-1596 USA\n"
      << "***\n"
      << "***  Permission is hereby granted, free of charge, to any person "
         "obtaining \n"
      << "***  a copy of this software and associated documentation files "
         "(the\n"
      << "***  Software, to deal in the Software without restriction, "
         "including\n"
      << "***  without limitation \n"
      << "***  the rights to use, copy, modify, merge, publish, distribute, "
         "sublicense, \n"
      << "***  and/or sell copies of the Software, and to permit persons to "
         "whom the \n"
      << "***  Software is furnished to do so, subject to the following "
         "conditions:\n"
      << "***\n"
      << "***  The above copyright notice and this permission notice shall be "
         "included\n"
      << "***  in all copies or substantial portions of the Software.\n"
      << "***\n"
      << "*** THE SOFTWARE IS PROVIDED AS IS, WITHOUT WARRANTY OF ANY KIND, \n"
      << "*** EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES\n"
      << "*** OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND "
         "NONINFRINGEMENT. \n"
      << "*** IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR "
         "ANY\n"
      << "*** CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF "
         "CONTRACT, TORT\n"
      << "*** OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE "
         "SOFTWARE OR\n"
      << "*** THE USE OR OTHER DEALINGS IN THE SOFTWARE.\n"
      << "***\n"
      << "***\n"
      << "*********************************************************************"
         "******\n"
      << endl;

  if (noParams.found() || helpRequest.found() || helpRequest1.found()) {
    cout << "  Use '-help' or '-f filename.aux' " << endl;
    cout << "  Saving options: \n"
         << "   -saveAsNodes     base_filename\n"
         << "   -saveAsCplace    base_filename\n"
         << "   -saveAsPlato     base_filename\n"
         << "   -saveLEFDEF      base_filename\n"
         << "   -markMacrosAsBlocks (for saveLEFDEF)\n"
         << "   -savePl  filename\n"
         << "  Misc options: \n"
         << "   -legal  (Use to remove any overlaps in existing placement)\n"
         << "   -snapToSites  (Use to snap cells to sites in existing "
            "placement)\n"
         << "   -spaceCellsEqually   (for each subrow, space cells equally)\n"
         << "  Plotting options: \n"
         << "   -xmin,-xmax,-ymin,-ymax  (for plotting range. Default means "
            "all)\n"
         << "   -plotNets        base_fileName\n"
         << "   -plotNodes       base_fileName\n"
         << "   -plotNodesWNames base_fileName\n"
         << "   -plot            base_fileName   (plots nodes and nets)\n"
         << "   -plotWNames      base_fileName   (the above + node names)\n"
         << "   -plotSites	  base_fileName   (plots the site map)\n"
         << "   -plotNodesWSites base_fileName   (the above + nodes)\n"
         << "   -plotRows        base_fileName   (plot the rows)\n"
         << endl;
    return 0;
  }
  cout << "RBPlacement Parameters" << endl;
  cout << rbParams << endl;

  StringParam auxFile("f", argc, argv);
  abkfatal(auxFile.found(), "must have -f <auxfilename>");
  RBPlacement rbplace(auxFile, rbParams);

  if (stdNodeHeights.found()) {
    double rowHeight = (rbplace.coreRowsBegin())->getHeight();
    const_cast<HGraphWDimensions &>(rbplace.getHGraph()).setNodeDims(rowHeight);
  }

  if (spaceCellsWCongInfo.found()) rbplace.spaceCellsWCongInfoInRows1();

  if (spaceCellsEqually.found()) rbplace.spaceCellsEquallyInRows();

  if (shiftCellsLeft.found()) rbplace.shiftCellsLeft();

  // add code that produces gnuplot files by traversing the hypergraph
  // and the layout from rbplace

  if (remSitesCong.found()) rbplace.removeSitesFromCongestedRgn(remSitesCong);

  MaxMem maxMem;

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

  if (remSitesMacro) rbplace.removeSitesBelowMacros();

  if (saveLEFDEF.found()) rbplace.saveLEFDEF(saveLEFDEF, markMacrosAsBlocks);

  if (saveAsNodes.found()) {
    if (!addNodesCong.found()) {
      if (fixMacros.found())
        rbplace.saveAsNodesNetsWts(saveAsNodes, fixMacros);
      else
        rbplace.saveAsNodesNetsWts(saveAsNodes);
    } else {
      rbplace.saveAsNodesNetsWtsWCongInfo(saveAsNodes, addNodesCong);
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

  if (saveAsNodesShredHW.found())
    rbplace.saveAsNodesNetsWtsShredHW(saveAsNodesShredHW);

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

  // cout<<endl<<endl<<auxFile<<" numMacros is :
  // "<<rbplace.getNumMacros()<<endl<<endl;

  return 0;
}
