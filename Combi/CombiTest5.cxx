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

// Gray-code partitioning traversal test.
// It steps through partitionings by Gray-code increments and prints the
// sequence for a requested size and number of parts.
// Use it when checking partitioning Gray-code generation.

#include <cstdlib>
#include <iostream>

#include "ABKCommon/abkassert.h"
#include "ABKCommon/abkio.h"
#include "grayPart.h"

using std::cout;

int main(int argc, char** argv) {
  abkfatal(argc >= 2, "Need size");
  unsigned size = atoi(argv[1]);
  unsigned numPart = 2;
  if (argc > 2) numPart = atoi(argv[2]);

  std::vector<unsigned> part(size, 0);
  cout << " Started with    : " << part;

  GrayCodeForPartitionings grayCode(size, numPart);

  while (!grayCode.finished()) {
    unsigned toIncr = grayCode.nextIncrement();
    cout << " Incrementing " << toIncr << "  : ";
    (++part[toIncr]) %= numPart;
    cout << part;
  }
}
