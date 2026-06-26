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

#ifndef _ORIENT_OPT_XFACE_H
#define _ORIENT_OPT_XFACE_H

#include "ABKCommon/abkcommon.h"
#include "RBPlace/RBPlacement.h"

class OrientOptParameters
{
    void setDefaults();
  public:
    Verbosity    verb;
    unsigned     numPassLimit;
    unsigned     cutOffGain;
    bool         minimizeCost;
    bool         checkSymmetry;
    bool         useSteiner;

    OrientOptParameters();
    OrientOptParameters(int argc, const char *argv[]);
    OrientOptParameters(const OrientOptParameters &srcParams);

    friend std::ostream& operator<<(std::ostream& os, const OrientOptParameters &params);
};


class OrientOptXFace;
typedef OrientOptXFace OrientOptimizerInterface;

class OrientOptXFace
{
  protected:

    OrientOptParameters   _params;
    std::vector<Orient>   _orients;
    unsigned              _nFlips;         // number of flipped cells
    double                _initialWL;
    double                _finalWL;

  public:

/*    OrientOptXFace(const DB& db, const OrientOptParameters& params); */
    OrientOptXFace(const RBPlacement& rbplace, const OrientOptParameters& params);

    const std::vector<Orient> &getOrients() const { return _orients; }
    unsigned getNumFlips() const { return _nFlips; }
    double getInitialWL()  const { return _initialWL; }
    double getFinalWL()    const { return _finalWL; }
};

#endif
