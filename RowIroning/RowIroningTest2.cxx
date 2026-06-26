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


#include "ABKCommon/abkcommon.h"
#include "ABKCommon/abkversion.h"
#include "DB/DB.h"
#include "RBPlFromDB/RBPlaceFromDB.h"
#include "NetTopology/hgWDimsFromDB.h"
#include "rowIroning.h"
#include "FixOr/greedyOrientOpt.h"
#include <sstream>
using std::stringstream;

void ri2d1pass(RBPlacement& rbplace, const RowIroningParameters& RIparams)
{
      RowIroningParameters TwoDRIparams;
      TwoDRIparams.numPasses = 1;
      TwoDRIparams.windowSize = 15;
      TwoDRIparams.overlap = 3;
      TwoDRIparams.twoDim = true;
      TwoDRIparams.mixed = false;
      TwoDRIparams.cluster = false;
      TwoDRIparams.useDP = true;
      TwoDRIparams.LRfirst = true;
      TwoDRIparams.verb = RIparams.verb;
      cout << TwoDRIparams;
      ironRows(rbplace, TwoDRIparams);
}

void ri2d2passes(RBPlacement& rbplace, const RowIroningParameters& RIparams)
{
      RowIroningParameters TwoDRIparams;
      TwoDRIparams.numPasses = 4;
      TwoDRIparams.windowSize = 10;
      TwoDRIparams.overlap = 3;
      TwoDRIparams.twoDim = true;
      TwoDRIparams.mixed = false;
      TwoDRIparams.cluster = false;
      TwoDRIparams.useDP = true;
      TwoDRIparams.LRfirst = true;
      TwoDRIparams.verb = RIparams.verb;
      cout << TwoDRIparams;
      ironRows(rbplace, TwoDRIparams);
}

void riBBpass(RBPlacement& rbplace, const RowIroningParameters& RIparams, unsigned windowSize)
{
      RowIroningParameters BBRIparams;
      BBRIparams.numPasses = 1;
      BBRIparams.windowSize = windowSize;
      BBRIparams.overlap = 2;
      BBRIparams.twoDim = false;
      BBRIparams.mixed = false;
      BBRIparams.cluster = false;
      BBRIparams.useDP = false;
      BBRIparams.LRfirst = true;
      BBRIparams.verb = RIparams.verb;
      cout << BBRIparams;
      ironRows(rbplace, BBRIparams);
}

void stdRowIroning(RBPlacement& rbplace,  const RowIroningParameters& RIparams)
{
      RowIroningParameters riparams(RIparams);
      riparams.windowSize = 16;
      cout << riparams;
      ironRows(rbplace, riparams);
}

void rowIroningPass(RBPlacement& rbplace, const RowIroningParameters& RIparams)
{
    vector<unsigned> rtype(5);
    for(unsigned loopCt = 0; loopCt < rtype.size(); ++loopCt)
        rtype[loopCt]=loopCt;

    RandomRawUnsigned r;
    random_shuffle(rtype.begin(), rtype.end(), r);      

    for(unsigned loopCt = 0; loopCt < rtype.size(); ++loopCt)
    {
        switch( rtype[loopCt] )
        {
            case 0:
                stdRowIroning(rbplace, RIparams);
                break;
            case 1:
                ri2d1pass(rbplace, RIparams);
                break;
            case 2:
                stdRowIroning(rbplace, RIparams);
                break;
            case 3:
                ri2d2passes(rbplace, RIparams);
                break;
            case 4:
                riBBpass(rbplace, RIparams, 8);
                break;
        }
    }
}


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

  

  RowIroningParameters RIparams(argc,argv);
  RandomRawUnsigned r;
  cout << "# The random seed for this run is " << RandomRoot::getExternalSeed() << endl;

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
      cout << " ====== Launching row ironing ... " << endl;
      cout << RIparams;
      
      cout <<"Initial Origin-to-Origin WL: "<<rbplace.evalHPWL(false)<<endl;
      double initialWL = rbplace.evalHPWL(true);
      cout <<"Initial Pin-to-Pin WL:       "<<initialWL<<endl;
     
      
      Timer RITimer;
      double wlImprovement = 1.0;
      double lastWL = rbplace.evalHPWL(true);
      double currWL = lastWL;
      unsigned numPasses = 0;
       
      while( wlImprovement > 0.05)
      {
        cout<<" ------ Beginning Composite Pass Number: " << numPasses << " ------ "<<endl;
        rowIroningPass(rbplace, RIparams);
        lastWL = currWL;
        currWL = rbplace.evalHPWL(true);
        
        wlImprovement = (-((currWL - lastWL)/lastWL))*100;

        //save after every composite pass
        stringstream ss;
        ss<<"out"<<static_cast<int>(numPasses)<<".pl";
        rbplace.savePlacement( ss.str().c_str() );
        rbplace.savePlacement( "out.riBest.pl" );
        numPasses++;

        //Stop after 24 hours guaranteed
        RITimer.split();
        double currTime = RITimer.getUserTime();
        if(currTime > 24 * 3600 ) wlImprovement = 0;
        RITimer.resume(); 
        
      }

      rbplace.remOverlaps();

      RITimer.stop();
      
      double finalWL = rbplace.evalHPWL(true);
      cout << " Final Origin-to-Origin WL:  "<<rbplace.evalHPWL(false)<<endl;
      cout << " Final Pin-to-Pin WL:        "<<finalWL<<endl;
      cout << " Total runtime:              "<<RITimer.getUserTime()<<endl;
      cout << " Overall improvement:        "<< (-((finalWL - initialWL)/initialWL)*100) <<"%"<<endl;
      cout << " Number of composite passes: "<< numPasses <<endl;
      
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



