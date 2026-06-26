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

#ifndef __LA_SPLITLARGE_H_
#define __LA_SPLITLARGE_H_

#include "Ctainers/bitBoard.h"
#include "baseBinSplitter.h"

class LookAheadSplitter : public BaseBinSplitter {
  unsigned _branchingFactor;
  unsigned _lookAheadLvls;

  std::vector<CapoBin*> _bTSLeftNeighbors;
  std::vector<CapoBin*> _bTSRightNeighbors;
  std::vector<CapoBin*> _bTSNeighborsAbove;
  std::vector<CapoBin*> _bTSNeighborsBelow;

  // bins created by the first partitioning...these are the
  // pairs of bins were' picking between
  std::vector<std::vector<CapoBin*> > _topLevelBins;
  std::vector<double> _estWL;
  std::vector<std::vector<unsigned> > _pinCounts;

  void attachNeighbors();

  BitBoard& _netIsInBin;

  void setNetsBitBoard();
  double evalBinsWL(const Placement& pl);

 public:
  LookAheadSplitter(CapoBin& binToSplit, CapoPlacer& capo,
                    unsigned branchingFactor, unsigned lookAheadLvls,
                    Verbosity verb = Verbosity("0_0_0"));

  virtual ~LookAheadSplitter();
};

#endif
