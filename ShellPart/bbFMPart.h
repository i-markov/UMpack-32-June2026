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

//////////////////////////////////////////////////////////////////////
//
// Created by Mike Oliver on Apr 24, 1999
// bbFMPart.h: interface for the BBFMPart class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(_BBFMPart_H_INCLUDED_)
#define _BBFMPart_H_INCLUDED_

#include "ABKCommon/abkassert.h"
#include "ABKCommon/paramproc.h"
#include "ABKCommon/verbosity.h"
#include "Partitioners/multiStartPart.h"
#include "SmallPart/bbPart.h"

class MaxMem;

class BBFMPart {
 public:
  class Parameters {
   public:
    Verbosity verb;
    unsigned fmStarts;
    BBPart::Parameters bbParams;
    MultiStartPartitioner::Parameters fmParams;

    Parameters() : verb("silent"), fmStarts(2), bbParams(), fmParams() {
      abkfatal(0, "disabling this ctor");
    }

    Parameters(int argc, const char *argv[])
        : verb(argc, argv),
          fmStarts(2),
          bbParams(argc, argv),
          fmParams(argc, argv) {
      UnsignedParam fms("fmStarts", argc, argv);
      if (fms.found()) fmStarts = fms;
    }

    Parameters(const Parameters &p)
        : verb(p.verb),
          fmStarts(p.fmStarts),
          bbParams(p.bbParams),
          fmParams(p.fmParams) {}
  };

 protected:
  PartitioningProblem &_problem;
  Parameters _params;
  bool _solvedExactly;
  BBPartBitBoardContainer &_bbBitBoards;

 public:
  BBFMPart(PartitioningProblem &problem, BBPartBitBoardContainer &bbBitBoards,
           MaxMem *maxMem, Parameters params = Parameters());
  virtual ~BBFMPart();

  bool solvedExactly() const { return _solvedExactly; }
};

#endif
