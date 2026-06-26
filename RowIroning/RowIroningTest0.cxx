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
#include "ABKCommon/abkmessagebuf.h"
#include "DB/DB.h"
#include "RBPlFromDB/RBPlaceFromDB.h"
#include "NetTopology/hgWDimsFromDB.h"
#include "rowIroning.h"
#include "FixOr/greedyOrientOpt.h"
#include "PlaceEvals/pinPlEval.h"

using std::cerr;
using std::cout;
using std::endl;


int main(int argc, const char *argv[])
{
  NoParams         noParams       (argc,argv);
  BoolParam        helpRequest    ("help",argc,argv);
  BoolParam        helpRequest1   ("h",   argc,argv);
  BoolParam        doFlipping     ("flip",argc,argv);
  StringParam      auxFileName    ("f",   argc,argv);

  cout << getABKMessageBuf();

  RowIroningParameters RIparams(argc,argv);
  RIparams.verb  = Verbosity("1 1 1");

  if (noParams.found() || helpRequest.found() || helpRequest1.found())
  {
     cout << " Use '-help' or '-f filename.aux' " << endl;
     return 0;
  }

  cout << RIparams;

  abkfatal(auxFileName.found(),"Usage: prog -f filename.aux <more params>"); 
  DB db(auxFileName);
  cerr<<"Design Name: "<<db.getDesignName()<<endl;

  RBPlaceParams rbParam(argc,argv);
  RBPlaceFromDB rbplace(db, rbParam);   

  cout <<"Initial Origin-to-Origin WL: "<<rbplace.evalHPWL(false)<<endl;
  cout <<"Initial Pin-to-Pin WL:       "<<rbplace.evalHPWL(true)<<endl;
#ifdef USEFLUTE
  bool usePins = true;
  if(RIparams.useSteiner)
    cout << "Routed WL estimate (FLUTE): "
         << pinPlEval::evalHPWL_Flute(rbplace.getPlacement(), rbplace.getHGraph(), usePins) << endl;
#endif

  SmallPlacementBitBoardContainer spBitBoards;
  MaxMem maxMem;

  if(!doFlipping.found())
  {
     ironRows(rbplace, RIparams, spBitBoards, &maxMem);
     rbplace.remOverlaps();
     cout << " Final Origin-to-Origin WL: "<<rbplace.evalHPWL(false)<<endl;
     cout << " Final Pin-to-Pin WL:       "<<rbplace.evalHPWL(true)<<endl;
#ifdef USEFLUTE
     if(RIparams.useSteiner)
       cout << "Routed WL estimate (FLUTE): "
            << pinPlEval::evalHPWL_Flute(rbplace.getPlacement(), rbplace.getHGraph(), usePins) << endl;
#endif
     db.setPlaceAndOrient(rbplace);
  }
  else
  {
     unsigned numPasses = RIparams.numPasses;
     RIparams.numPasses = 1;
     RIparams.verb = "0_0_0";

     for(unsigned p = 0; p < numPasses; p++)
     {
	 cout<<" Ironing/Flipping Pass: "<<p<<endl;

         Timer ironTime;
         ironRows(rbplace, RIparams, spBitBoards, &maxMem);
         ironTime.stop();
         cout <<"   RI Runtime:  " << ironTime.getUserTime()<<endl;
         cout <<"   RI Origin-to-Origin WL: "<<rbplace.evalHPWL(false)<<endl;
         cout <<"   RI Pin-to-Pin WL:       "<<rbplace.evalHPWL(true)<<endl;
#ifdef USEFLUTE
         if(RIparams.useSteiner)
           cout << "Routed WL estimate (FLUTE): "
                << pinPlEval::evalHPWL_Flute(rbplace.getPlacement(), rbplace.getHGraph(), usePins) << endl;
#endif

         OrientOptParameters flipParams;
         Timer flipTimer;
         GreedyOrientOpt flipper(rbplace,flipParams);
         flipTimer.stop();

         rbplace<<flipper.getOrients();

         db.setPlaceAndOrient(rbplace);

         cout <<"   FLIP Runtime: "<<flipTimer.getUserTime()
		<<" (flipped "<<flipper.getNumFlips()<<" cells )"<<endl;
         cout <<"   FLIP Origin-to-Origin WL: "<<rbplace.evalHPWL(false)<<endl;
         cout <<"   FLIP Pin-to-Pin WL:       "<<rbplace.evalHPWL(true)<<endl;
#ifdef USEFLUTE
         if(RIparams.useSteiner)
           cout << "Routed WL estimate (FLUTE): "
                << pinPlEval::evalHPWL_Flute(rbplace.getPlacement(), rbplace.getHGraph(), usePins) << endl;
#endif
     }

     cout << endl;
     cout <<"Final Origin-to-Origin WL: "<<rbplace.evalHPWL(false)<<endl;
     cout <<"Final Pin-to-Pin WL:       "<<rbplace.evalHPWL(true)<<endl;
#ifdef USEFLUTE
     if(RIparams.useSteiner)
       cout << "Routed WL estimate (FLUTE): "
            << pinPlEval::evalHPWL_Flute(rbplace.getPlacement(), rbplace.getHGraph(), usePins) << endl;
#endif
  }


  db.saveDEF("newDef");
  return 0;
}
