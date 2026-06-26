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

/*
  This subHGraph constructor is for use with BoxPlacer.  It creates
  a HGraph from:
   ->vector<CapoBin*> 	all CapoBins
   ->HGraphFixed		hgraphOfNetlist

 The hgraph is a copy of the netlist hgraph, with the addition of:
  4*#CapoBins nodes
  4*#cells in all CapoBins edges.

  for each CapoBin there are 4 new nodes, and an edge from each cell in
  the bin to each of those nodes.
*/

#ifndef __CAPO_CONSTINGHG_H__
#define __CAPO_CONSTINGHG_H__

#include "HGraph/subHGraph.h"
#include "capoBin.h"

class ConstrainingHGraph : public HGraphFixed {
  friend class CapoPlacer;

  const std::vector<CapoBin*>& _allBins;
  const HGraphFixed& _netlistHG;

 public:
  ConstrainingHGraph(const std::vector<CapoBin*>& allBins,
                     const HGraphFixed& hgraphOfNetlist,
                     const Verbosity verb = Verbosity("1_1_1"));

  virtual ~ConstrainingHGraph() {}
};

#endif
