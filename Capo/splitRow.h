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

// SplitRowOfBins takes a vector of bins (a contig row w same y-span),
// combines then into one bin, splits this horizontally, and then
// constructs 2*orig# bins new bins that have as equal as possible
// whitespace distrubition

#ifndef __SPLITROW_H_
#define __SPLITROW_H_

#include "baseBinSplitter.h"

class SplitRowOfBins : public BaseBinSplitter {
  std::vector<CapoBin*>& _rowOfBins;
  static bit_vector isInTopBin;

  bool _wasGoodSplit;

  void buildSubBins(std::vector<CapoBin*>& bigBins);

  void buildOneSubRow(CapoBin* subRowBin, bool isTopBin,
                      std::vector<CapoBin*>& newBins);

 public:
  SplitRowOfBins(std::vector<CapoBin*>& rowOfBins, CapoPlacer& capo,
                 CapoSplitterParams params = CapoSplitterParams());

  virtual ~SplitRowOfBins() {}

  bool wasGoodSplit() const { return _wasGoodSplit; }
};

#endif
