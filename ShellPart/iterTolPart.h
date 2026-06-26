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

// Created by Andrew Caldwell on Apr 27, 1999

// iterative partitioner calls FM several times, with decreasing tolerance
//  if given a problem where there are movable cells if size > 1/2 the
//  original tolerance.

#ifndef __SHELLPART_ITERTOLPART_H__
#define __SHELLPART_ITERTOLPART_H__

#include "ABKCommon/paramproc.h"
#include "ABKCommon/verbosity.h"
#include "FMPart/fmPart.h"

#include <vector>

class MaxMem;

class IterativeTolPart {
 public:
  class Parameters {
   public:
    Verbosity verb;

    double initTolRatio;  // these control the tolerance given
                          // to the first partitioning instance.
    double initTolPerct;  // the tol is set to
                          //  max(origTol, initTolRatio*largestCellArea,
                          //  initTolPerct)

    unsigned firstFMNumPasses;  // number of FM passes to do at the
                                // the first FM call (if there are
                                // going to be two of them).
                                // A value of 0 means a complete FM run.

    bool firstFMUseClip;  // if true, clip will be used on the
                          // FIRST FM call, only

    MultiStartPartitioner::Parameters fmParams;

    Parameters()
        : verb("silent"),
          initTolRatio(3.0),
          initTolPerct(20),
          firstFMNumPasses(0),
          firstFMUseClip(false),
          fmParams() {}

    Parameters(int argc, const char* argv[])
        : verb(argc, argv),
          initTolRatio(3.0),
          initTolPerct(20),
          firstFMNumPasses(0),
          firstFMUseClip(false),
          fmParams(argc, argv) {
      DoubleParam initTRat("initTolRatio", argc, argv);
      DoubleParam initTPct("initTolPerct", argc, argv);
      UnsignedParam numP("firstFMNumPasses", argc, argv);
      BoolParam FMUseClip("firstFMUseClip", argc, argv);

      if (initTRat.found()) initTolRatio = initTRat;
      if (initTPct.found()) initTolPerct = initTPct;
      if (numP.found()) firstFMNumPasses = numP;
      firstFMUseClip = FMUseClip.found();
    }

    friend std::ostream& operator<<(std::ostream& os, const Parameters& params);
  };

 protected:
  PartitioningProblem& _problem;
  Parameters _params;
  std::vector<std::vector<double> > _origMaxCap;
  std::vector<std::vector<double> > _origMinCap;

 public:
  IterativeTolPart(PartitioningProblem& problem, MaxMem* maxMem,
                   Parameters params);
};

#endif
