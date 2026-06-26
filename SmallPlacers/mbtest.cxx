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

// Small helper sanity test for the `mbitset` utility.
// It exercises construction, setting, unsetting, and indexing so the data
// structure can be checked quickly from the command line.

#include "mbset.h"

int main() {
  mbitset b1(0, 0), b2(1, 0), b3(0, 1), b4(2, 2), b5(3, 7);
  cout << b1 << b2 << b3 << b4 << b5;
  mbitset s1(0, 0);
  cout << "set bit 3" << endl;
  s1.set(3);
  cout << s1;
  cout << "set bit 32" << endl;
  s1.set(32);
  cout << s1;
  cout << "unset bit 2" << endl;
  s1.unset(2);
  cout << s1;
  cout << "unset bit 3" << endl;
  s1.unset(3);
  cout << s1;
  cout << "test bit 3 : " << s1[3] << endl;
  cout << "test bit 32: " << s1[32] << endl;

  return 0;
}
