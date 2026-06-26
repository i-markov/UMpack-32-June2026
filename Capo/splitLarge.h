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

// created by Andrew Caldwell on 10/17/99 caldwell@cs.ucla.edu

// SplitLarge contains classes that are designed to partiton
// small CapoBins.  Large is expected to mean > 200 nodes.

// The following 2 classes are defined:
// SplitLargeBinHorizontally(CapoBin&, const MLParams&, Verbosity)
// SplitLargeBinVertically  (CapoBin&, const MLParams&, Verbosity)

// each take as input a reference to a CapoBin, and populate
// vector of pointers to new CapoBins.

#ifndef __SPLITLARGE_H_
#define __SPLITLARGE_H_

#include <vector>

#include "baseBinSplitter.h"

enum SplitPtFixedType { SPLITPT_UNFIXED, SPLITPT_FIXED, SPLITPT_HINT };

class SplitLargeBinHorizontally : public BaseBinSplitter {
  friend class PartitioningProblemForCapo;

  unsigned _numMLSets;

  void callPartitioner(PartitioningProblemForCapo& problem,
                       PlacementWOrient* placement = NULL);
  bool callMLPart(PartitioningProblemForCapo& problem,
                  PlacementWOrient* placement = NULL);

 public:
  SplitLargeBinHorizontally(CapoBin& binToSplit, CapoPlacer& capo,
                            CapoSplitterParams params = CapoSplitterParams(),
                            PlacementWOrient* placement = NULL,
                            SplitPtFixedType splitPtFixedType = SPLITPT_UNFIXED,
                            double ySplit = 0.);

  virtual ~SplitLargeBinHorizontally() {}
};

class SplitLargeBinVertically : public BaseBinSplitter {
  friend class PartitioningProblemForCapo;

  unsigned _numMLSets;

  void callPartitioner(PartitioningProblemForCapo& problem,
                       PlacementWOrient* placement = NULL);
  bool callMLPart(PartitioningProblemForCapo& problem,
                  PlacementWOrient* placement = NULL);

 public:
  SplitLargeBinVertically(CapoBin& binToSplit, CapoPlacer& capo,
                          CapoSplitterParams params = CapoSplitterParams(),
                          PlacementWOrient* placement = NULL,
                          SplitPtFixedType splitPtFixedType = SPLITPT_UNFIXED,
                          double xSplit = 0.);

  virtual ~SplitLargeBinVertically() {}
};

void string_vec_to_char_star_view(const std::vector<std::string>& in,
                                  std::vector<const char*>& out);

#endif
