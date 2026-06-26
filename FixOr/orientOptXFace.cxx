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


// Created 990210 by Stefanus Mantik, VLSI CAD ABKGROUP, UCLA/CS

#include "orientOptXFace.h"

using std::ostream;
using std::endl;

void OrientOptParameters::setDefaults()
{
    numPassLimit = UINT_MAX;
    cutOffGain = 0;
    minimizeCost = true;
    checkSymmetry = true;
    useSteiner = false;
}

OrientOptParameters::OrientOptParameters()
    : verb("1 1 1")
{
    setDefaults();
}

OrientOptParameters::OrientOptParameters(int argc, const char* argv[])
{
    setDefaults();

    StringParam    verbIn("verb", argc, argv);
    UnsignedParam  numPassLimitIn("nPass", argc, argv);
    UnsignedParam  cutOffGainIn("cutOff", argc, argv);
    BoolParam      maximizeCost("maximize", argc, argv);

    if(verbIn.found()) verb = Verbosity(verbIn);
    if(numPassLimitIn.found())  numPassLimit = numPassLimitIn;
    if(cutOffGainIn.found()) cutOffGain = cutOffGainIn;
    if(maximizeCost) minimizeCost = false;
}

OrientOptParameters::OrientOptParameters(const OrientOptParameters &srcParams)
    : verb(srcParams.verb),
      numPassLimit(srcParams.numPassLimit),
      cutOffGain(srcParams.cutOffGain),
      minimizeCost(srcParams.minimizeCost),
      checkSymmetry(srcParams.checkSymmetry),
      useSteiner(srcParams.useSteiner)
{}

ostream& operator<<(ostream& os, const OrientOptParameters &params)
{
    os << "Number of pass limit = " << params.numPassLimit << endl
       << "Cut off gain = " << params.cutOffGain << endl;
    if(!params.minimizeCost) os << "Use Maximize cost" << endl;
    if(params.checkSymmetry) os << "Check Symmetry" << endl;
    else os << "Do Not Check Symmetry" << endl;
    return os;
}

/*

OrientOptXFace::OrientOptXFace(const DB& db, const OrientOptParameters& params) :
      _params(params), _orients(db.spatial.orientations), _nFlips(0), _initialWL(0), _finalWL(0)
{}

*/

OrientOptXFace::OrientOptXFace(const RBPlacement& rbplace, const OrientOptParameters& params) :
      _params(params), _orients(rbplace.getPlacement().getSize(),Orient("N")), _nFlips(0), _initialWL(0), _finalWL(0)
{ 
  const Orient *orients = rbplace.getPlacement().getOrients();
  for(unsigned i = 0; i < _orients.size(); i++)
  {
     _orients[i] = orients[i];
  }
}
