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
***  EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
***  MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
***  IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
***  CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
***  TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
***  SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
***
***************************************************************************/

// Correlation/clumping regression driver.
// It compares the effect of clumping and unclumping around the center on
// placement correlation and rank correlation metrics.

// Created by Igor Markov, September 1997, VLSICAD ABKGroup UCLA/CS

#include <fstream>
#include <iostream>
using std::cout;
using std::endl;
#include "plStats.h"
#include "transf.h"

int main() {
  // Placement pl(Placement::_Grid1,3);
  Placement pl(Placement::_Random, 100);
  BBoxToBBox(BBox(0, 0, 1, 1), BBox(2, 2, 10, 10))(pl);
  // ofstream("orig.gen.pla") << pl;
  {
    Placement pl1(pl);
    ClumpInCenter clump(pl1.getBBox(), Point(5, 5), 3.0);  // alpha==3.0
    clump(pl1);
    cout << "Clumping points in pl near the center -- producing pl1" << endl;
    cout << "PlCorr(pl,pl1)      = " << PlCorr(pl, pl1) << endl;
    cout << "PlRankCorr(pl,pl1)  = " << PlRankCorr(pl, pl1) << endl;
    //    ofstream("clumped.pla") << pl1;
  }
  {
    Placement pl2(pl);
    ClumpInCenter unclump(pl2.getBBox(), Point(5, 5),
                          1.0 / 3.0);  // alpha==1/3.0
    unclump(pl2);
    cout << "Unclumping points in pl from the center -- producing pl2" << endl;
    cout << "PlCorr(pl,pl2)      = " << PlCorr(pl, pl2) << endl;
    cout << "PlRankCorr(pl,pl2)  = " << PlRankCorr(pl, pl2) << endl;
    //    ofstream("unclumped.pla") << pl2;
  }
}
