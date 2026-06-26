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

// Created by Igor Markov, Sat Jul 21, 2001

#ifndef __HG_DELAY_INFO_H_
#define __HG_DELAY_INFO_H_

#include <iostream>
#include <vector>

class HGDelayInfo {
 public:
  bool populated;
  // The Cint_per_unit and Rint_per_unit are functions of netlength:
  //   Cint_per_unit = exp(Cint_K * log(length) + Cint_C)
  //   Rint_per_unit = exp(Rint_K * log(length) + Rint_C)

  double rInt_c, rInt_k;
  double cInt_c, cInt_k;
  // these two arrays use pin indexing from HGraphWDims
  std::vector<double> drivRes;  // per pin, only meaningful for driver pins
  std::vector<double> sinkCap;  // per pin, only meaningful for sink   pins

  HGDelayInfo(unsigned numPins = 0)
      : populated(false),
        rInt_c(0),
        rInt_k(0),
        cInt_c(0),
        cInt_k(0),
        drivRes(numPins, 0.0),
        sinkCap(numPins, 0.0) {}
};

std::ostream& operator<<(std::ostream& out, const HGDelayInfo& di);

#endif
