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

// TrivialStats and FreqDistr walkthrough.
// It demonstrates the expected-minimums helper on two samples and prints a
// summary-statistics report for the second sample.
// Use it as a compact example of the basic stats utilities.

#include <algorithm>
#include <iostream>

#include "ABKCommon/paramproc.h"
#include "linRegre.h"
#include "stats.h"

using std::cin;
using std::copy;
using std::cout;
using std::endl;

int main(int argc, char *argv[]) {
  UnsignedParam num("n", argc, argv);
  /*
      std::vector<double> x,y;
      double x_in, y_in;
      for(unsigned i = 0; i < num; i++)
      {
          cin >> x_in >> y_in;
          x.push_back(x_in);
          y.push_back(y_in);
      }

      LinearRegression lreg(x,y);
      cout << "LINREG: C = " << lreg.getC() << ", K = " << lreg.getK() << endl;
      Correlation corr(x,y);
      RankCorr rcorr(x,y);
      cout << "Correlation: " << corr << endl;
      cout << "Rank Correlation: " << rcorr << endl;
      TrivialStatsWithStdDev tx(x), ty(y);
      cout << tx << "\nSTD=" << tx.getStdDev()
           <<"\n===\n" << ty << "\nSTD=" << ty.getStdDev() << endl;

      std::vector<double> err;
      for(i = 0; i < x.size(); i++)
      {
          err.push_back(abs(y[i]-x[i]*lreg.getK()-lreg.getC()));
      }
      TrivialStatsWithStdDev ers(err);
      cout << ers;
  */
  std::vector<double> x;
  double x_in;
  for (unsigned i = 0; i < num; i++) {
    cin >> x_in;
    x.push_back(x_in);
  }
  TrivialStatsWithStdDev tx(x);
  cout << tx;

  return 0;
}
