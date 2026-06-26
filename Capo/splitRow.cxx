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

#ifdef _MSC_VER
#pragma warning(disable : 4786)
#endif

#include "splitRow.h"

#include <algorithm>

#include "splitLarge.h"

using std::cout;
using std::endl;
using std::sort;
using std::string;
using std::vector;

bit_vector SplitRowOfBins::isInTopBin;

SplitRowOfBins::SplitRowOfBins(vector<CapoBin*>& rowOfBins, CapoPlacer& capo,
                               CapoSplitterParams params)
    : BaseBinSplitter(*rowOfBins[0], capo, params),
      _rowOfBins(rowOfBins),
      _wasGoodSplit(true) {
  abkfatal(_rowOfBins.size() > 1, "SplitRowOfBins expects > 1 bin");

  if (_params.verb.getForMajStats() > 2)
    cout << "Using RowSplitter on a row of " << rowOfBins.size() << " bins"
         << endl;

  unsigned k;
  for (k = 0; k < rowOfBins.size(); k++) {
    if (_params.verb.getForMajStats() > 3)
      cout << *rowOfBins[k] << endl;
    else if (_params.verb.getForMajStats() > 2)
      cout << "[" << rowOfBins[k]->getIndex() << "]" << endl;
  }

  if (isInTopBin.size() < capo.getNetlistHGraph().getNumNodes())
    isInTopBin = bit_vector(capo.getNetlistHGraph().getNumNodes());

  // construct the big bin
  sort(rowOfBins.begin(), rowOfBins.end(), CompareSrcBinsByX());
  CapoBin rowBin(rowOfBins);
  const_cast<CapoPlacer&>(_capo).updateInfoAboutBin(&rowBin);

  if (_params.verb.getForMajStats() > 4)
    cout << "  Created RowBin \n" << rowBin << endl << endl;

  SplitLargeBinHorizontally rowBinSplitter(rowBin, capo, _params);

  buildSubBins(rowBinSplitter.getNewBins());

  if (_wasGoodSplit) {
    for (k = 0; k < rowOfBins.size(); k++) rowOfBins[k]->unLinkNeighbors();
    for (k = 0; k < _newBins.size(); k++) _newBins[k]->linkNeighbors();
  } else  // put it back the way it was
  {
    for (k = 0; k < rowOfBins.size(); k++)
      const_cast<CapoPlacer&>(_capo).updateInfoAboutBin(rowOfBins[k]);
    for (k = 0; k < _newBins.size(); k++) {
      delete _newBins[k];
      _newBins[k] = 0;
    }
    _newBins.clear();
  }
}

void SplitRowOfBins::buildSubBins(vector<CapoBin*>& bigBins) {
  CapoBin* topBin = bigBins[0];
  CapoBin* botBin = bigBins[1];

  vector<unsigned>::const_iterator cId;
  for (cId = topBin->cellIdsBegin(); cId != topBin->cellIdsEnd(); cId++)
    isInTopBin[*cId] = true;
  for (cId = botBin->cellIdsBegin(); cId != botBin->cellIdsEnd(); cId++)
    isInTopBin[*cId] = false;

  vector<CapoBin*> newTopBins;
  buildOneSubRow(topBin, true, newTopBins);

  vector<CapoBin*> newBotBins;
  buildOneSubRow(botBin, false, newBotBins);

  // hook up the neighbors, sibblings and names

  _newBins = newTopBins;
  _newBins.insert(_newBins.end(), newBotBins.begin(), newBotBins.end());
}

void SplitRowOfBins::buildOneSubRow(CapoBin* subRowBin, bool isTopBin,
                                    vector<CapoBin*>& newBins) {
  double partLayoutArea = subRowBin->getCapacity();

  // layout area to allocate per unit of cell area
  double layoutAreaRatio = partLayoutArea / subRowBin->getTotalCellArea();

  vector<CapoBin*> noNeighbor;  // for top and bottom neighbors..currently
                                // unpopulated!!
  vector<CapoBin*> lNeighbors;
  vector<CapoBin*> rNeighbors;

  vector<unsigned> cellsInBin;
  double areaInBin;
  double leftEdge = subRowBin->getBBox().xMin;
  double rightEdge;

  vector<double> newBinCaps;

  // total area of all bins left of the current one.
  double cumLayoutArea = 0;

  const HGraphWDimensions& hgraph = _capo.getNetlistHGraph();
  unsigned b;
  for (b = 0; b < _rowOfBins.size(); b++) {
    const CapoBin& oldBin = *_rowOfBins[b];
    cellsInBin.clear();
    areaInBin = 0;

    vector<unsigned>::const_iterator cId;
    for (cId = oldBin.cellIdsBegin(); cId != oldBin.cellIdsEnd(); cId++)
      if ((isInTopBin[*cId] && isTopBin) || (!isInTopBin[*cId] && !isTopBin)) {
        cellsInBin.push_back(*cId);
        areaInBin += hgraph.getWeight(*cId);
      }

    if (areaInBin < oldBin.getTotalCellArea() * 0.2) {
      cout << "Part Area: " << areaInBin
           << " total bin area: " << oldBin.getTotalCellArea() << endl;
      _wasGoodSplit = false;
    }
    if (areaInBin == 0) continue;  // don't produce 0 area bins

    cumLayoutArea += (areaInBin * layoutAreaRatio);

    if (b < _rowOfBins.size() - 1)  // not the end...find the xSplit
    {
      double p0WSWeight = 1.;
      double p1WSWeight = 1.;
      rightEdge =
          subRowBin->findXSplit(cumLayoutArea, partLayoutArea - cumLayoutArea,
                                0, p0WSWeight, p1WSWeight, newBinCaps);
    } else
      rightEdge = subRowBin->getBBox().xMax;

    BBox partBoundary = subRowBin->getBBox();
    partBoundary.xMin = leftEdge;
    partBoundary.xMax = rightEdge;

    string name(oldBin.getName());
    if (isTopBin)
      name += "R0";
    else
      name += "R1";

    CapoBin* newBin =
        new CapoBin(cellsInBin, subRowBin->getRows(), partBoundary, noNeighbor,
                    noNeighbor, lNeighbors, rNeighbors, _capo, name);

    if (_params.verb.getForActions() > 4)
      cout << "\nBuilt new bin: \n" << *newBin << endl;

    leftEdge = rightEdge;
    newBins.push_back(newBin);
  }
}
