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

// BitBoard smoke test.
// It sets, clears, and queries bits while printing the internal bit-vector
// representation.
// Use it to verify the basic bit-board container behavior.

#include <iostream>

#include "bitBoard.h"

using std::cout;
using std::endl;

int main(int argc, char** argv) {
  BitBoard bb(10);
  const bit_vector& bv = bb.getBoardSpace();
  cout << bb;
  unsigned k;
  for (k = 0; k != bv.size(); k++) cout << bv[k];
  cout << endl;
  cout << " Bit 5 is set to " << bb.isBitSet(5) << endl;
  cout << " Set bits 3 and 7 " << endl;
  bb.setBit(3);
  bb.setBit(7);
  cout << bb;
  for (k = 0; k != bv.size(); k++) cout << bv[k];
  cout << endl;
  cout << " Set bit 3 once again " << endl;
  bb.setBit(3);
  cout << bb;
  cout << " Clear the BitBoard " << endl;
  bb.clear();
  cout << bb;
  for (k = 0; k != bv.size(); k++) cout << bv[k];
  cout << endl;
}
