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

// Interval-sequence smoke test.
// It exercises complement, merge, and intersection on a few small interval
// sets.
// Use it to verify the core IntervalSeq operations.

#include "interval.h"

using std::cout;
using std::endl;

int main() {
  IntervalSeq iSeq1;

  iSeq1.push_back(Interval(-10, 10));
  iSeq1.push_back(Interval(20, 30));

  cout << "IntervalSequence 1: " << endl << iSeq1 << endl;
  iSeq1.complement();
  cout << "Complemented:       " << endl << iSeq1 << endl;
  iSeq1.complement();
  cout << "Original:           " << endl << iSeq1 << endl;

  IntervalSeq iSeq2;
  iSeq2.push_back(Interval(-5, 5));
  iSeq2.push_back(Interval(8, 22));
  iSeq2.push_back(Interval(25, 40));

  cout << "IntervalSequence 2: " << endl << iSeq2 << endl;

  IntervalSeq tempISeq1 = iSeq1;

  tempISeq1.merge(iSeq2);
  cout << "ISeq1 Union ISeq2:     " << endl << tempISeq1 << endl;
  iSeq1.intersect(iSeq2);
  cout << "ISeq1 Intersect ISeq2: " << endl << iSeq1 << endl;
}
