/**************************************************************************
***    
*** Copyright (c) 2026 IonQ, Inc.
*** Copyright (c) 1995-2000 Regents of the University of California,
***               Andrew E. Caldwell, Andrew B. Kahng and Igor L. Markov
*** Copyright (c) 2000-2012 Regents of the University of Michigan,
***               Saurabh N. Adya, Jarrod A. Roy, David A. Papa and
***               Igor L. Markov
***
***  Contact author(s): abk@cs.ucsd.edu, imarkov@umich.edu
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


// RowIroning regression driver (Bookshelf aux-file path).
//
// Reads a row-based placement from -f <design.aux>, optionally removes
// overlaps (skip with -skipLegal), then calls ironRows() to slide windows
// of cells within (and optionally across) placement rows and accept moves
// that reduce half-perimeter wirelength.  Reports initial and final
// pin-to-pin HPWL.  Key options:
//   -smPlAlgo Branch | DynamicP   small-placer algorithm (default Branch)
//   -ironWindow N                  window size (default 12)
//   -ironPasses N                  number of passes (default 1)
//   -ironTwoDim                    enable two-row ironing
//   -noRowIroning                  skip ironing (legalization only)
//   -save                          write final placement to out.pl
//   -plot* <base>                  gnuplot visualization scripts
//
// This test exercises the pure-RBPlacement code path; it does not require
// the DB or RBPlFromDB packages.

#include "ABKCommon/abkcommon.h"
#include "ABKCommon/abkversion.h"
#include "ABKCommon/abkmessagebuf.h"
#include "rowIroning.h"
#include "FixOr/greedyOrientOpt.h"
#include "PlaceEvals/pinPlEval.h"
using std::cout;
using std::endl;


int main(int argc, const char *argv[])
{
  NoParams         noParams       (argc,argv);
  BoolParam        helpRequest    ("help",argc,argv);
  BoolParam        helpRequest1   ("h",   argc,argv);
  BoolParam        doFlipping     ("flip",argc,argv);
  StringParam      auxFileName    ("f",   argc,argv);
  BoolParam        save           ("save",argc,argv);
  BoolParam        noRowIroning   ("noRowIroning",argc,argv);

  StringParam   plotNets("plotNets", argc, argv);
  StringParam   plotNodes("plotNodes", argc, argv);
  StringParam   plotNodesAndNets("plot", argc, argv);
  StringParam   plotNodesWNames("plotNodesWNames", argc, argv);
  StringParam   plotNodesAndNetsWNames("plotWNames", argc, argv);
  StringParam   plotSites("plotSites", argc, argv);
  StringParam   plotNodesWSites("plotNodesWSites", argc, argv);  
  BoolParam     keepOverlaps    ("skipLegal",argc,argv);

  cout << getABKMessageBuf();

  RowIroningParameters RIparams(argc,argv);

  if (noParams.found() || helpRequest.found() || helpRequest1.found())
  {
     cout << " Use '-help' or '-f filename.aux' " << endl;
     return 0;
  }

  abkfatal(auxFileName.found(),"Usage: prog -f filename.aux <more params>"); 

  RBPlaceParams rbParam(argc,argv);
  RBPlacement rbplace(auxFileName, rbParam);   

  if(!keepOverlaps.found())
    {
      rbplace.remOverlaps();
    }

  if(!noRowIroning.found())
    {
      SmallPlacementBitBoardContainer spBitBoards;
      cout << " ====== Launching row ironing ... " << endl;
      cout << RIparams;
      
      cout <<"Initial Origin-to-Origin WL: "<<rbplace.evalHPWL(false)<<endl;
      cout <<"Initial Pin-to-Pin WL:       "<<rbplace.evalHPWL(true)<<endl;
#ifdef USEFLUTE
      bool usePins = true;
      if(RIparams.useSteiner)
        cout << " Routed WL estimate (FLUTE): "
             << pinPlEval::evalHPWL_Flute(rbplace.getPlacement(), rbplace.getHGraph(), usePins) << endl;
#endif

      MaxMem maxMem;      
      Timer RITimer;

      ironRows(rbplace, RIparams, spBitBoards, &maxMem);
      rbplace.remOverlaps();

      RITimer.stop();
      
      cout << " Final Origin-to-Origin WL: "<<rbplace.evalHPWL(false)<<endl;
      cout << " Final Pin-to-Pin WL:       "<<rbplace.evalHPWL(true)<<endl;
#ifdef USEFLUTE
      if(RIparams.useSteiner)
        cout << "Routed WL estimate (FLUTE): "
             << pinPlEval::evalHPWL_Flute(rbplace.getPlacement(), rbplace.getHGraph(), usePins) << endl;
#endif
      
    }

  Plotters::RBPlacePlotter::Parameters allPlotParams(argc,argv,false);
  if(plotNets.found()) {
    Plotters::RBPlacePlotter::Parameters plotParams(allPlotParams);
    plotParams.plotNets=true;
    rbplace.saveAsPlot(plotParams,plotNets);
  }

  if(plotNodes.found()) {
    Plotters::RBPlacePlotter::Parameters plotParams(allPlotParams);
    plotParams.plotNodes=true;
    rbplace.saveAsPlot(plotParams,plotNodes);
  }

  if(plotNodesAndNets.found()) {
    Plotters::RBPlacePlotter::Parameters plotParams(allPlotParams);
    plotParams.plotNodes=true;
    plotParams.plotNets=true;
    rbplace.saveAsPlot(plotParams,plotNodesAndNets);
  }

  if(plotNodesWNames.found()) {
    Plotters::RBPlacePlotter::Parameters plotParams(allPlotParams);
    plotParams.plotNodes=true;
    plotParams.plotNodesNames=true;
    rbplace.saveAsPlot(plotParams,plotNodesWNames);
  }

  if(plotNodesAndNetsWNames.found()) {
    Plotters::RBPlacePlotter::Parameters plotParams(allPlotParams);
    plotParams.plotNodes=true;
    plotParams.plotNodesNames=true;
    plotParams.plotNets=true;
    rbplace.saveAsPlot(plotParams,plotNodesAndNetsWNames);
  }

  if(plotSites.found()) {
    Plotters::RBPlacePlotter::Parameters plotParams(allPlotParams);
    plotParams.plotSiteMap=true;
    rbplace.saveAsPlot(plotParams,plotSites);
  }

  if(plotNodesWSites.found()) {
    Plotters::RBPlacePlotter::Parameters plotParams(allPlotParams);
    plotParams.plotSiteMap=true;
    plotParams.plotNodes=true;
    rbplace.saveAsPlot(plotParams,plotNodesWSites);
  }

  if(save.found())
    {
      cout << "Saving out.pl" << endl;
      rbplace.savePlacement("out.pl");
    }
  
  return 0;
}



