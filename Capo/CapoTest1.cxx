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

// Regression test driver for Capo.
// It exercises representative package behavior with sample inputs or
// command-line options and prints results for the regression harness.
// Individual files vary the data set, but follow the same compare-against-
// baseline pattern.

// Flexible Bookshelf placement driver for Capo.
//
// Normal mode reads the design from -f <file.aux>, runs Capo, and can write:
//   -savePl <file.pl>          final Bookshelf placement
//   -plot* <baseFileName>      gnuplot visualization scripts
//
// Ordering mode is selected by -saveOrder <file>.  It reads only the netlist,
// forces unit-size nodes, builds a synthetic single-row placement with one
// slot per movable node, runs Capo, and writes movable node names in final
// left-to-right row order.  -savePl and -plot* can be used in the same run to
// inspect the resulting one-dimensional placement.
#ifdef _MSC_VER
#pragma warning(disable : 4786)
#endif

#include <fstream>
#include <iostream>
#include <vector>

#include "ABKCommon/abkassert.h"
#include "ABKCommon/abkmessagebuf.h"
#include "ABKCommon/abktimer.h"
#include "ABKCommon/abkversion.h"
#include "ABKCommon/verbosity.h"
#include "RBPlace/RBPlacement.h"
#include "capoPlacer.h"

using std::cerr;
using std::cout;
using std::endl;
using std::flush;
using std::ofstream;
using std::vector;

static RBPlacement *buildOneRowPlacement(const char *auxFileName,
                                         const RBPlaceParams &rbParam) {
  // Build the synthetic layout directly so unit dimensions are in place before
  // row capacity is computed.
  RBPlaceParams netlistParams(rbParam);
  netlistParams.makeAllSrcSnk = false;
  HGraphWDimensions *hgraph =
      new HGraphWDimensions(auxFileName, NULL, NULL, netlistParams);

  const unsigned numNodes = hgraph->getNumNodes();
  const unsigned numTerms = hgraph->getNumTerminals();
  const unsigned numCoreCells = numNodes - numTerms;
  abkfatal(numCoreCells != 0, "No core cells to place");

  for (unsigned i = 0; i < numNodes; ++i) {
    hgraph->setNodeWidth(1.0, i);
    hgraph->setNodeHeight(1.0, i);
    hgraph->updateNodeMacroInfo(i, false);
  }

  PlacementWOrient placement(numNodes, Point(-1.0, -1.0));
  Orient nOrient("N");
  for (unsigned i = 0; i < numTerms; ++i) {
    const bool onLeft = (i % 2 == 0);
    placement[i] = Point(onLeft ? -1.0 : static_cast<double>(numCoreCells), 0);
    placement.setOrient(i, nOrient);
  }
  for (unsigned i = numTerms; i < numNodes; ++i) {
    placement[i] = Point(static_cast<double>(i - numTerms), 0.0);
    placement.setOrient(i, nOrient);
  }

  RBPSite site(1.0, 1.0, Symmetry("Y"));
  vector<RBPCoreRow> coreRows;
  coreRows.push_back(RBPCoreRow(0.0, nOrient, site, placement, 1.0));
  coreRows.back().appendNewSubRow(0.0, static_cast<double>(numCoreCells));

  return new RBPlacement(hgraph, placement, coreRows, rbParam);
}

static void saveOrdering(const RBPlacement &rbplace, const char *fileName) {
  // The populated subrow cell vectors represent final row order.
  ofstream out(fileName);
  abkfatal(out, "Could not open ordering output file");

  const HGraphWDimensions &hgraph = rbplace.getHGraph();
  unsigned numWritten = 0;

  for (itRBPCoreRow row = rbplace.coreRowsBegin();
       row != rbplace.coreRowsEnd(); ++row) {
    for (itRBPSubRow subRow = row->subRowsBegin();
         subRow != row->subRowsEnd(); ++subRow) {
      for (itRBPCellIds cell = subRow->cellIdsBegin();
           cell != subRow->cellIdsEnd(); ++cell) {
        if (*cell < hgraph.getNumTerminals()) continue;
        out << hgraph.getNodeNameByIndex(*cell) << endl;
        ++numWritten;
      }
    }
  }

  abkfatal(numWritten == hgraph.getNumNodes() - hgraph.getNumTerminals(),
           "Ordering output does not contain every core cell");
}

int main(int argc, const char *argv[]) {
  NoParams noParams(argc, argv);
  BoolParam helpRequest("help", argc, argv);
  BoolParam helpRequest1("h", argc, argv);
  StringParam savePl("savePl", argc, argv);
  StringParam saveOrder("saveOrder", argc, argv);
  StringParam plotNodes("plotNodes", argc, argv);
  StringParam plotNodesAndNets("plot", argc, argv);
  StringParam plotNodesWNames("plotNodesWNames", argc, argv);
  StringParam plotNodesAndNetsWNames("plotWNames", argc, argv);
  StringParam plotSites("plotSites", argc, argv);
  StringParam plotNodesWSites("plotNodesWSites", argc, argv);
  StringParam plotRows("plotRows", argc, argv);
  // BoolParam        doCellFlipping ("flip",argc,argv);
  // BoolParam        doCellSpacing  ("cellSpace",argc,argv);

  Verbosity verb(argc, argv);

  cout << "Capo 11.0 (June 2026)" << endl;
  cout << "Compiled " << __DATE__ << " at " << __TIME__ << endl << endl;

  cout
      << "\n(c) Copyright, 1998-2000 Regents of the University of California\n";
  cout << "(c) Copyright, 2000-2007 Regents of the University of Michigan\n\n";
  cout << " Authors: Saurabh N. Adya, Andrew E. Caldwell, Andrew B. Kahng\n";
  cout << "          Igor L. Markov, David A. Papa and Jarrod A. Roy\n";
  cout << " Email:   imarkov@eecs.umich.edu\n";
  cout << " http://vlsicad.eecs.umich.edu/BK/PDtools/" << endl << endl;

  cout << getABKMessageBuf() << flush;

  PRINT_VERSION_INFO

  cout << endl
       << "Permission is hereby granted, without written agreement and \n"
       << "without license or royalty fee, to use, copy, modify, and \n"
       << "distribute and sell this software and its documentation for \n"
       << "any purpose, provided that the above copyright notice, this \n"
       << "permission notice, and the following two paragraphs appear \n"
       << "in all copies of this software as well as in all copies of \n"
       << "supporting documentation.\n\n";
  cout << "THIS SOFTWARE AND SUPPORTING DOCUMENTATION ARE PROVIDED \"AS IS\".\n"
       << "The Microelectronics Advanced Research Corporation (MARCO), the\n"
       << "Gigascale Silicon Research Center (GSRC),  and \n"
       << "(\"PROVIDERS\") MAKE NO WARRANTIES, whether express \n"
       << "or implied, including warranties of merchantability or fitness\n"
       << "for a particular purpose or noninfringement, with respect to \n"
       << "this software and supporting documentation.\n"
       << "Providers have NO obligation to provide ANY support, assistance,\n"
       << "installation, training or other services, updates, enhancements\n"
       << "or modifications related to this software and supporting\n"
       << "documentation.\n\n";
  cout << "Providers shall NOT be liable for ANY costs of procurement of \n"
       << "substitutes, loss of profits, interruption of business, or any \n"
       << "other direct, indirect, special, consequential or incidental \n"
       << "damages arising from the use of this software and its \n"
       << "documentation, whether or not Providers have been advised of \n"
       << "the possibility of such damages." << endl
       << endl;

  cout << CmdLine(argc, argv);

  // ------------- Process parameters before getting bogged down in computations

  bool saveMem = false;
  CapoParameters capoParams(argc, argv, saveMem);
  MaxMem maxMem;
  capoParams.maxMem = &maxMem;

  if (noParams.found() || helpRequest.found() || helpRequest1.found()) {
    cout << "  Use '-help' or '-f filename.aux' " << endl;
    cout << "  To save output, use -savePl filename " << endl;
    cout << "  To place a netlist in one row and save a node ordering, use"
         << " -saveOrder filename" << endl;
    cout << "  To save a gnuplot script, use one of:" << endl;
    cout << "    -plotNodes <baseFileName>" << endl;
    cout << "    -plotNodesWNames <baseFileName>" << endl;
    cout << "    -plot <baseFileName>" << endl;
    cout << "    -plotWNames <baseFileName>" << endl;
    cout << "    -plotSites <baseFileName>" << endl;
    cout << "    -plotNodesWSites <baseFileName>" << endl;
    cout << "    -plotRows <baseFileName>" << endl;
    return 0;
  }

  cout << capoParams << endl;

  StringParam auxFileName("f", argc, argv);
  abkfatal(auxFileName.found(), "Usage: prog -f filename.aux <more params>");

  RBPlaceParams rbParam(argc, argv);
  cout << rbParam << endl;
  RBPlacement *rbplace = saveOrder.found()
                             ? buildOneRowPlacement(auxFileName, rbParam)
                             : new RBPlacement(auxFileName, rbParam);

  Timer capoTimer;
  CapoPlacer capo(*rbplace, capoParams);
  capoTimer.stop();

  cout << "After Capo  " << endl;
  cout << "  Capo Runtime: " << capoTimer.getUserTime() << endl;

  cout << "  RBPl Est HPWL: " << rbplace->evalHPWL(false) << endl;
  ;
  cout << "  RBPl Est Pin-to-Pin HPWL: " << rbplace->evalHPWL(true) << endl;
  ;
  cout << "  RBPl Est  WWL: " << rbplace->evalWeightedWL(false) << endl;
  ;
  cout << "  RBPl Est Pin-to-Pin  WWL: " << rbplace->evalWeightedWL(true)
       << endl;
  ;

  if (savePl.found()) rbplace->savePlacement(savePl);
  if (saveOrder.found()) saveOrdering(*rbplace, saveOrder);

  Plotters::RBPlacePlotter::Parameters allPlotParams(argc, argv, false);
  if (plotNodes.found()) {
    Plotters::RBPlacePlotter::Parameters plotParams(allPlotParams);
    plotParams.plotNodes = true;
    rbplace->saveAsPlot(plotParams, plotNodes);
  }
  if (plotNodesAndNets.found()) {
    Plotters::RBPlacePlotter::Parameters plotParams(allPlotParams);
    plotParams.plotNodes = true;
    plotParams.plotNets = true;
    rbplace->saveAsPlot(plotParams, plotNodesAndNets);
  }
  if (plotNodesWNames.found()) {
    Plotters::RBPlacePlotter::Parameters plotParams(allPlotParams);
    plotParams.plotNodes = true;
    plotParams.plotNodesNames = true;
    rbplace->saveAsPlot(plotParams, plotNodesWNames);
  }
  if (plotNodesAndNetsWNames.found()) {
    Plotters::RBPlacePlotter::Parameters plotParams(allPlotParams);
    plotParams.plotNodes = true;
    plotParams.plotNodesNames = true;
    plotParams.plotNets = true;
    rbplace->saveAsPlot(plotParams, plotNodesAndNetsWNames);
  }
  if (plotSites.found()) {
    Plotters::RBPlacePlotter::Parameters plotParams(allPlotParams);
    plotParams.plotSiteMap = true;
    rbplace->saveAsPlot(plotParams, plotSites);
  }
  if (plotNodesWSites.found()) {
    Plotters::RBPlacePlotter::Parameters plotParams(allPlotParams);
    plotParams.plotSiteMap = true;
    plotParams.plotNodes = true;
    rbplace->saveAsPlot(plotParams, plotNodesWSites);
  }
  if (plotRows.found()) {
    Plotters::RBPlacePlotter::Parameters plotParams(allPlotParams);
    plotParams.plotRows = true;
    rbplace->saveAsPlot(plotParams, plotRows);
  }

  delete rbplace;
  return 0;
}
