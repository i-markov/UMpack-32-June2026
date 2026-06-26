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

// Summary-statistics smoke test.
// It prints TrivialStats and TrivialStatsWithStdDev for a few small vectors.
// Use it to verify the basic descriptive-statistics helpers.

#include <iostream>

#include "ABKCommon/abkconst.h"
#include "stats.h"

using std::cout;
using std::endl;

int main(int argc, char *argv[]) {
  (void)argc;
  (void)argv;
  std::vector<double> dataSet1(3, Pi);
  std::vector<unsigned> dataSet2(3, 8);

  std::vector<double> dataSet3(4);
  std::vector<unsigned> dataSet4(4);

  // std::iota(dataSet3.begin(),dataSet3.end(),0);
  {
    for (unsigned i = 0; i < dataSet3.size(); i++) dataSet3[i] = i;
  }

  // std::iota(dataSet4.begin(),dataSet4.end(),0);
  {
    for (unsigned i = 0; i < dataSet4.size(); i++) dataSet4[i] = i;
  }

  cout << " Stats1  " << TrivialStats(dataSet1) << " Stats1a "
       << TrivialStatsWithStdDev(dataSet1) << " Stats2  "
       << TrivialStats(dataSet2) << " Stats2a "
       << TrivialStatsWithStdDev(dataSet2) << " Stats3  "
       << TrivialStats(dataSet3) << " Stats3a "
       << TrivialStatsWithStdDev(dataSet3) << " Stats4  "
       << TrivialStats(dataSet4) << " Stats4a "
       << TrivialStatsWithStdDev(dataSet4);

  return 0;
}
