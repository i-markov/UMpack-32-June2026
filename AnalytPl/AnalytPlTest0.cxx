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

// Analytical placement regression driver.
// It loads an RBPlacement from -f, builds the analytical solver over the
// movable cells, and runs a placement update so the regression harness can
// compare the resulting coordinates and optional plots.

#include "ABKCommon/abkassert.h"
#include "ABKCommon/abktimer.h"
#include "ABKCommon/paramproc.h"
#include "analytPl.h"

using std::cout;
using std::endl;
using std::max;
using std::vector;

int main(int argc, const char* argv[]) {
  NoParams noParams(argc, argv);
  BoolParam helpRequest("help", argc, argv);
  BoolParam helpRequest1("h", argc, argv);
  StringParam auxFile("f", argc, argv);
  UnsignedParam skipNetsLargerThan("skipNetsLargerThan", argc, argv);
  BoolParam usePinOffsets("usePinOffsets", argc, argv);
  BoolParam noCOG("noCOG", argc, argv);
  BoolParam useLinear("useLinear", argc, argv);
  UnsignedParam numLinearIter("numLinearIter", argc, argv);

  StringParam plotNets("plotNets", argc, argv);
  StringParam plotNodes("plotNodes", argc, argv);
  StringParam plotNodesWNames("plotNodesWNames", argc, argv);
  StringParam plotNodesAndNets("plot", argc, argv);

  if (noParams.found() || helpRequest.found() || helpRequest1.found()) {
    cout << "  Use '-help' or '-f filename.aux' " << endl;
    cout << "Options\n"
         << "-usePinOffsets <bool>\n"
         << "-skipNetsLargerThan <unsigned>\n"
         << "-noCOG <bool>\n"
         << "-useLinear <bool>\n"
         << "-numLinearIter <unsigned> [5]\n\n"
         << "-plotNodes       base_fileName\n"
         << "-plotNodesWNames base_fileName\n"
         << "-plot            base_fileName   (plots nodes and nets)\n\n";
    exit(0);
  }

  RBPlacement::Parameters rbParams(argc, argv);

  abkfatal(auxFile.found(), "must have -f <auxfilename>");

  RBPlacement rbplace(auxFile, rbParams);
  BBox layoutBBox = rbplace.getBBox(true);
  const HGraphWDimensions& hgWDims = rbplace.getHGraph();

  // layoutBBox.ShrinkTo(0.8);

  vector<vector<unsigned> > nodeIds(1);
  vector<BBox> binBBox(1);
  for (unsigned i = 0; i < hgWDims.getNumNodes(); ++i) {
    if (!hgWDims.isTerminal(i)) nodeIds[0].push_back(i);
  }
  binBBox[0].xMin = layoutBBox.xMin;
  binBBox[0].xMax = layoutBBox.xMax;
  binBBox[0].yMin = layoutBBox.yMin;
  binBBox[0].yMax = layoutBBox.yMax;

  // double change=DBL_MAX;
  // double binMaxChange=0;

  AnalyticalSolver solver(rbplace, nodeIds, binBBox);
  solver.initNodeLocsToBinCenter();
  if (skipNetsLargerThan.found())
    solver.setSkipNetsLargerThan(skipNetsLargerThan);
  if (usePinOffsets.found()) solver.setUsePinOffsets(true);

  unsigned numLinIter = 5;
  if (numLinearIter.found()) numLinIter = numLinearIter;

  double maxDimension = max(layoutBBox.getHeight(), layoutBBox.getWidth());
  double epsilon = maxDimension / 100;

  Timer totalTime;
  if (!noCOG.found()) {
    if (useLinear.found())
      solver.solveLinearMin(epsilon, numLinIter);
    else
      solver.solveQuadraticMin(epsilon);
  } else {
    if (useLinear.found())
      solver.solveLinearSOR(epsilon, 0, numLinIter);
    else
      solver.solveSOR(epsilon, 0);
  }

  totalTime.stop();
  cout << "Time taken = " << totalTime.getUserTime() << " seconds\n";

  vector<vector<Point> >& nodeLocs = solver.getNodeLocs();
  for (unsigned binIdx = 0; binIdx < nodeIds.size(); ++binIdx) {
    for (unsigned i = 0; i < nodeIds[binIdx].size(); ++i) {
      rbplace.updatePlacement(nodeIds[binIdx][i], nodeLocs[binIdx][i]);
    }
  }

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

  if (plotNodesWNames.found()) {
    Plotters::RBPlacePlotter::Parameters plotParams(allPlotParams);
    plotParams.plotNodes = true;
    plotParams.plotNodesNames = true;
    rbplace.saveAsPlot(plotParams, plotNodesWNames);
  }

  if (plotNodesAndNets.found()) {
    Plotters::RBPlacePlotter::Parameters plotParams(allPlotParams);
    plotParams.plotNodes = true;
    plotParams.plotNets = true;
    rbplace.saveAsPlot(plotParams, plotNodesAndNets);
  }

  return 0;
}
