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

// Prim MST benchmark with file input.
// It reads a point set from `-in`, builds the MST, and prints the cost and
// elapsed time.
// Use it for simple performance comparisons on point files.

// Created by Igor Markov, July 1997, VLSICAD ABKGroup UCLA/CS

#include <vector>

#include "ABKCommon/abkio.h"
#include "ABKCommon/abktimer.h"
#include "ABKCommon/paramproc.h"
#include "primMST.h"

using std::cout;
using std::endl;
using std::flush;
using std::ifstream;
using std::vector;

int main(int argc, char* argv[]) {
  StringParam inputFileName("in", argc, argv);
  vector<Point> points;
  ifstream input(inputFileName);
  int lineNo = 0;
  double x, y;
  while (!input.eof()) {
    input >> eathash(lineNo);
    if (input.eof()) break;
    input >> x >> y >> skiptoeol;
    points.push_back(Point(x, y));
  }

  Timer mstTime;
  PrimMST mst(points);
  mstTime.stop();
  cout << "MST Cost : " << mst.getCost() << endl
       << "Time: " << mstTime << endl
       << flush;

  return 0;
}
