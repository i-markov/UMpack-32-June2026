/**************************************************************************
***    
*** Copyright (c) 2026 IonQ, Inc.
*** Copyright (c) 1995-2000 Regents of the University of California,
***               Andrew E. Caldwell, Andrew B. Kahng and Igor L. Markov
*** Copyright (c) 2000-2012 Regents of the University of Michigan,
***               Saurabh N. Adya, Jarrod A. Roy, David A. Papa and
***               Igor L. Markov
***
***  Contact author(s): abk@cs.ucsd.edu, imarkov@umich.edu
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


// Created 990210 by Stefanus Mantik, VLSI CAD ABKGROUP, UCLA/CS

#include "greedyOrientOpt.h"
#include "heapGainCtr.h"
#include "Ctainers/bitBoard.h"
#include "ABKCommon/abkcommon.h"
#ifdef USEFLUTE
#include "Flute/flute.h"
#endif
#include <map>

using std::map;
using std::min;
using std::max;
using std::pair;
using std::cout;
using std::endl;
using std::vector;

/*

GreedyOrientOpt::GreedyOrientOpt(const DB& db,
                                 const OrientOptParameters& params) :
    OrientOptimizerInterface(db,params)
{
    unsigned numNets = db.getNetlist().getNumNets(),
             numCells = db.getNetlist().getNumCells();
    unsigned i;
    vector<BBox>   netBBox(numNets), minPossibleBBox(numNets);
    vector<double> oldCost(numCells, 0);
    vector<double> newCost(numCells, 0);
    vector<double> gain(numCells, -DBL_MAX);
    vector<vector<unsigned> > needUpdate(numCells);
    BitBoard       flipableCells(numCells);
    vector<Orient>& ors = _orients;
    vector<pair<Point, Point> > pinLocations(db.getNetlist().getNumPins());

    itCellGlobal c;
    for(c = db.getNetlist().cellsBegin(); c != db.getNetlist().cellsEnd(); c++)
    {
        if(!(*c)->isCore()) {
          continue;    // only consider core cells
        }
        unsigned symm = static_cast<unsigned>((*c)->getMaster().getSymmetry());
        if(!(symm == 2 || symm == 3 || symm == 6 || symm == 7)) {
            continue;        // ignore the cell if it is not Y symmetry
        }
        flipableCells.setBit((*c)->getIndex());
    }
    if(flipableCells.getNumBitsSet() == 0) 
    {
        abkwarn(0, "no orientation change is permitted");
        return;
    }

    // Store the actual location of pins for flipped and unflipped
    itPinLocal p;
    for(c = db.getNetlist().cellsBegin();
        c != db.getNetlist().cellsEnd(); c++)
    {
        unsigned cellIdx = (*c)->getIndex();
        for(p = (*c)->pinsBegin(); p != (*c)->pinsEnd(); p++)
            pinLocations[(*p)->getIndex()].first = db.locatePin(**p, ors);
        if(flipableCells.isBitSet(cellIdx))
        {
            ors[cellIdx].flipHoriz();
            for(p = (*c)->pinsBegin(); p != (*c)->pinsEnd(); p++)
                pinLocations[(*p)->getIndex()].second = db.locatePin(**p, ors);
            ors[cellIdx].flipHoriz();
        }
        else
            for(p = (*c)->pinsBegin(); p != (*c)->pinsEnd(); p++)
                pinLocations[(*p)->getIndex()].second = 
                    pinLocations[(*p)->getIndex()].first;
    }

    // Store net BBox and minimum possible BBox
    itNetGlobal n;
    for(n = db.getNetlist().netsBegin(); n != db.getNetlist().netsEnd(); n++)
    {
        unsigned netIdx = (*n)->getIndex();
        BBox nBox;
        double minX = DBL_MAX, maxX = -DBL_MAX;
        for(p = (*n)->pinsBegin(); p != (*n)->pinsEnd(); p++)
        {
            const Point& fPt = pinLocations[(*p)->getIndex()].first,
                       & sPt = pinLocations[(*p)->getIndex()].second;
            nBox += fPt; //pinLocations[(*p)->getIndex()].first;
            if(fPt.x < sPt.x)
            {
                minX  = min(minX, sPt.x);
                maxX  = max(maxX, fPt.x);
            }
            else
            {
                minX  = min(minX, fPt.x);
                maxX  = max(maxX, sPt.x);
            }
        }
        netBBox[netIdx] = nBox;
        nBox.xMin = minX;
        nBox.xMax = maxX;
        minPossibleBBox[netIdx] = nBox;
    }

    // Build needUpdate list
    for(n = db.getNetlist().netsBegin(); n != db.getNetlist().netsEnd(); n++)
    {
        unsigned netIdx = (*n)->getIndex();
        BBox nBox = minPossibleBBox[netIdx];
        double halfPerim = netBBox[netIdx].getHalfPerim();
        _initialWL += halfPerim;
                   // looking for dependency between nets end cells
                   // this is for checking whether when a cell flipped,
                   // an update to the net is necessary or not
        itCellLocal lc;
        for(lc = (*n)->cellsBegin(); lc != (*n)->cellsEnd(); lc++)
        {
            unsigned cellIdx = (*lc)->getIndex();
            abkfatal(cellIdx < oldCost.size(), "Index too much");
            if(!flipableCells.isBitSet(cellIdx)) continue;
            oldCost[cellIdx] += halfPerim;
            Point rightEdge, leftEdge;
            rightEdge = leftEdge = db.spatial.locations[cellIdx];
            rightEdge.x += (*lc)->getWidth();   // only X axis 
                                               // because we only flip on Y
            if(nBox.contains(leftEdge) && nBox.contains(rightEdge))
                continue;
            abkfatal(cellIdx < needUpdate.size(), "Index too much");
            needUpdate[cellIdx].push_back(netIdx);
        }
    }

    for(c = db.getNetlist().cellsBegin();
        c != db.getNetlist().cellsEnd(); c++)
    {
        unsigned cellIdx = (*c)->getIndex();
        if (!flipableCells.isBitSet(cellIdx)) continue;
        newCost[cellIdx] = oldCost[cellIdx];
        for(i = 0; i < needUpdate[cellIdx].size(); i++)
        {
            unsigned netIdx = needUpdate[cellIdx][i];
            newCost[cellIdx] -= netBBox[netIdx].getHalfPerim();
            const dbNet& net = db.getNetlist().getNetByIdx(netIdx);
            BBox nBox;
            for(p = net.pinsBegin(); p != net.pinsEnd(); p++)
                if((*p)->getCell().getIndex() == cellIdx)
                    nBox += pinLocations[(*p)->getIndex()].second;
                else
                    nBox += pinLocations[(*p)->getIndex()].first;
            newCost[cellIdx] += nBox.getHalfPerim();
        }
        if(_params.minimizeCost)
            gain[cellIdx] = oldCost[cellIdx] - newCost[cellIdx];
        else
            gain[cellIdx] = newCost[cellIdx] - oldCost[cellIdx];
    }
    BitBoard visitedCells(numCells);
    HeapGainContainer  heapGain(gain, flipableCells);

                   // Optimization loop starts here
    while(gain[heapGain.getTop()] > _params.cutOffGain)
    {
        visitedCells.clear();
        unsigned cellIdx = heapGain.getTop();
        ors[cellIdx].flipHoriz();
        _nFlips++;

                          // problem with swap so use "old-style" one
        double swapBucket = oldCost[cellIdx];
        oldCost[cellIdx] = newCost[cellIdx];
        newCost[cellIdx] = swapBucket;
        gain[cellIdx] = -gain[cellIdx];
        heapGain.heapify(cellIdx);

        if(_params.verb.getForActions() > 1)
            cout << "delta : " << oldCost[cellIdx] - newCost[cellIdx] << endl;
        for(i = 0; i < needUpdate[cellIdx].size(); i++)
        {
            unsigned netIdx = needUpdate[cellIdx][i];
            const dbNet& net = db.getNetlist().getNetByIdx(netIdx);
            BBox nBox;
            for(p = net.pinsBegin(); p != net.pinsEnd(); p++)
                nBox += db.locatePin(**p, ors);
            double delta = nBox.getHalfPerim() -
                           netBBox[netIdx].getHalfPerim();
            netBBox[netIdx] = nBox;
            itCellLocal lc;
            for(lc = net.cellsBegin(); lc != net.cellsEnd(); lc++)
            {
                unsigned lCellIdx = (*lc)->getIndex();
                if(!flipableCells.isBitSet(lCellIdx)) continue;
                if(lCellIdx == cellIdx) continue;
                oldCost[lCellIdx] += delta;
                if(visitedCells.isBitSet(lCellIdx)) continue;
                visitedCells.setBit(lCellIdx);
            }
        }
        const vector<unsigned>& 
            updateCellList = visitedCells.getIndicesOfSetBits();

        for(i = 0; i < updateCellList.size(); i++)
        {
            unsigned lCellIdx = updateCellList[i];
            ors[lCellIdx].flipHoriz();
            newCost[lCellIdx] = oldCost[lCellIdx];
            for(unsigned j = 0; j < needUpdate[lCellIdx].size(); ++j)
            {
                unsigned netIdx = needUpdate[lCellIdx][j];
                const dbNet& net = db.getNetlist().getNetByIdx(netIdx);
                BBox nBox;
                for(p = net.pinsBegin(); p != net.pinsEnd(); p++)
                    nBox += db.locatePin(**p, ors);
                newCost[lCellIdx] += nBox.getHalfPerim() -
                                     netBBox[netIdx].getHalfPerim();
            }
            ors[lCellIdx].flipHoriz();
            if(_params.minimizeCost)
                gain[lCellIdx] = oldCost[lCellIdx] - newCost[lCellIdx];
            else
                gain[lCellIdx] = newCost[lCellIdx] - oldCost[lCellIdx];
            heapGain.heapify(lCellIdx);
        }
    }
    for(i = 0; i < numNets; i++)
        _finalWL += netBBox[i].getHalfPerim();
}

*/

GreedyOrientOpt::GreedyOrientOpt(const RBPlacement& rbplace,
                                 const OrientOptParameters& params) :
    OrientOptimizerInterface(rbplace,params)
{
    unsigned numNets = rbplace.getHGraph().getNumEdges(),
             numCells = rbplace.getHGraph().getNumNodes();
    vector<BBox> netBBox(numNets);
    vector<double> netCosts(numNets,0.);
    vector<Orient>& ors = _orients;
    vector< vector<unsigned> > edge_id_to_pin(numNets);
    BitBoard flipableCells(numCells);

#ifdef USEFLUTE
    if(_params.useSteiner)
    {
      readLUT();
    }
#endif

    itHGFEdgeGlobal e;
    itHGFNodeLocal c,n;
    unsigned nodeoffset;
    _initialWL = 0.;
    for(e = rbplace.getHGraph().edgesBegin(); e != rbplace.getHGraph().edgesEnd(); ++e)
    {
      const unsigned &netIdx = (*e)->getIndex();
      edge_id_to_pin[netIdx] = vector<unsigned>((*e)->getDegree(),0);
      BBox nBox;
      for(nodeoffset = 0, n = (*e)->nodesBegin(); n != (*e)->nodesEnd(); ++nodeoffset, ++n)
      {
        const unsigned &cellIdx = (*n)->getIndex();
        BBox pinBBox = rbplace.getHGraph().nodesOnEdgePinOffsetBBox(nodeoffset,netIdx,ors[cellIdx]);
        pinBBox.TranslateBy(rbplace[cellIdx]);
        nBox.expandToInclude(pinBBox);
      }
      netBBox[netIdx] = nBox;
      netCosts[netIdx] = getNetCost(rbplace, netIdx, nBox);
      if(!nBox.isEmpty())
      {
        _initialWL += nBox.getHalfPerim();
      }
    }

    unsigned edgeoffset, pin_id = 0;
    for(c = rbplace.getHGraph().nodesBegin(); c != rbplace.getHGraph().nodesEnd(); c++)
    {
        const unsigned &cellIdx = (*c)->getIndex();
        for(edgeoffset = 0, e = (*c)->edgesBegin(); e != (*c)->edgesEnd(); ++edgeoffset, ++e)
        {
            for(nodeoffset = 0, n = (*e)->nodesBegin();
                n != (*e)->nodesEnd() && (*n)->getIndex() != cellIdx; ++nodeoffset, ++n);
            abkfatal(n != (*e)->nodesEnd(),"could not find corresponding node from edge");
            edge_id_to_pin[(*e)->getIndex()][nodeoffset] = pin_id;
            pin_id++;
        }
    }

    // call Vert if macros
    if(rbplace.getNumMacros() > 0)
    {
      optimizeVert(rbplace,netBBox,netCosts,flipableCells,edge_id_to_pin);
      flipableCells.clear();
    }

    // call Horiz
    optimizeHoriz(rbplace,netBBox,netCosts,flipableCells,edge_id_to_pin);

    _finalWL = 0.;
    for(unsigned i = 0; i < numNets; i++)
    {
      if(!netBBox[i].isEmpty())
      {
        _finalWL += netBBox[i].getHalfPerim();
      }
    }
}

typedef pair<double,double> Span;

#define steinerCutoff 20

double GreedyOrientOpt::getNetCost(const RBPlacement& rbplace, unsigned netId, const BBox &netBBox) const
{
   const HGFEdge& net = rbplace.getHGraph().getEdgeByIdx(netId);
   unsigned netDegree = net.getDegree();
   if(netDegree <= 1)
   {
     return 0.;
   }
   else if(netDegree <= 3)
   {
     return netBBox.getHalfPerim();
   }
#ifdef USEFLUTE
   else if(_params.useSteiner && netDegree <= steinerCutoff)
   {
     vector<Point> points;
     unsigned nodeoffset;
     itHGFNodeLocal node;
     for(nodeoffset = 0, node = net.nodesBegin(); node != net.nodesEnd(); ++nodeoffset, ++node)
     {
       BBox pinOffset = rbplace.getHGraph().nodesOnEdgePinOffsetBBox(nodeoffset,netId,_orients[(*node)->getIndex()]);
       pinOffset.TranslateBy(rbplace[(*node)->getIndex()]);
       points.push_back(pinOffset.getBottomLeft());
       if(!equalDouble(pinOffset.xMin,pinOffset.xMax) || !equalDouble(pinOffset.yMin,pinOffset.yMax))
       {
         points.push_back(pinOffset.getTopRight());
       }
     }

     std::sort(points.begin(), points.end());
     vector<Point>::iterator new_end = std::unique(points.begin(),points.end());
     points.erase(new_end,points.end());

     double *x = new double[points.size()];
     double *y = new double[points.size()];
     
     for(unsigned i = 0; i < points.size(); ++i)
     {
       x[i] = points[i].x;
       y[i] = points[i].y;
     }

     Tree fltree = flute(points.size(), x, y, ACCURACY);
     double length = fltree.length;

     free(fltree.branch);

     delete [] x;
     delete [] y;

     return length;
   }
#endif
   else
   {
     return netBBox.getHalfPerim();
   }
}

void GreedyOrientOpt::optimizeVert(const RBPlacement& rbplace, vector<BBox>& netBBox, vector<double> &netCosts,
BitBoard& flipableCells, const vector< vector<unsigned> >& edge_id_to_pin)
{
    unsigned numNets = rbplace.getHGraph().getNumEdges(),
             numCells = rbplace.getHGraph().getNumNodes();
    vector<double> oldCost(numCells, 0.);
    vector<double> newCost(numCells, 0.);
    vector<double> gain(numCells, -DBL_MAX);
    vector<Span> minPossibleYSpan(numNets,Span(DBL_MAX,-DBL_MAX));
    vector<Span> maxYSpans(numCells,Span(DBL_MAX,-DBL_MAX));
    vector< vector<unsigned> > needUpdate(numCells);
    vector<Orient>& ors = _orients;
    vector< pair<Span, Span> > pinLocations(rbplace.getHGraph().getNumPins());

    itHGFNodeGlobal c;
    for(c = rbplace.getHGraph().nodesBegin(); c != rbplace.getHGraph().nodesEnd(); c++)
    {
        if(!rbplace.isCoreCell((*c)->getIndex())) {
          continue;    // only consider core cells
        }
        if(rbplace.isFixed((*c)->getIndex())) {
          continue;
        }
        if(!rbplace.getHGraph().isNodeMacro((*c)->getIndex())) {
          continue;
        }
        Symmetry symm = rbplace.getHGraph().getNodeSymmetry((*c)->getIndex());
        if(_params.checkSymmetry && !symm.x) {
            continue;        // ignore the cell if it is not X symmetry
        }
        flipableCells.setBit((*c)->getIndex());
    }
    if(flipableCells.getNumBitsSet() == 0)
    {
        abkwarn(0, "no vert orientation change is permitted");
        return;
    }

    // Store the actual location of pins for flipped and unflipped
    unsigned edgeoffset, nodeoffset;
    itHGFEdgeLocal e;
    itHGFNodeLocal node;
    unsigned pin_id = 0;
    for(c = rbplace.getHGraph().nodesBegin(); c != rbplace.getHGraph().nodesEnd(); ++c)
    {
        unsigned cellIdx = (*c)->getIndex();
        for(edgeoffset = 0, e = (*c)->edgesBegin(); e != (*c)->edgesEnd(); ++edgeoffset, ++e)
        {
            BBox pinOffset = rbplace.getHGraph().edgesOnNodePinOffsetBBox(edgeoffset,cellIdx,ors[cellIdx]);
            pinOffset.TranslateBy(rbplace[cellIdx]);
            pinLocations[pin_id].first.first  = pinOffset.yMin;
            pinLocations[pin_id].first.second = pinOffset.yMax;
            if(flipableCells.isBitSet(cellIdx))
            {
              ors[cellIdx].flipVert();
              pinOffset = rbplace.getHGraph().edgesOnNodePinOffsetBBox(edgeoffset,cellIdx,ors[cellIdx]);
              pinOffset.TranslateBy(rbplace[cellIdx]);
              pinLocations[pin_id].second.first  = pinOffset.yMin;
              pinLocations[pin_id].second.second = pinOffset.yMax;
              ors[cellIdx].flipVert();
            }
            else
            {
              pinLocations[pin_id].second = pinLocations[pin_id].first;
            }
            ++pin_id;
        }
    }

    // Store minimum and maximum possible YSpans
    itHGFEdgeGlobal n;
    for(n = rbplace.getHGraph().edgesBegin(); n != rbplace.getHGraph().edgesEnd(); ++n)
    {
        unsigned netIdx = (*n)->getIndex();
        double minY = DBL_MAX, maxY = -DBL_MAX;
        for(nodeoffset = 0, node = (*n)->nodesBegin(); node != (*n)->nodesEnd(); ++nodeoffset, ++node)
        {
            //unsigned pinId = edge_id_to_pin.find(pair<unsigned,unsigned>(netIdx,nodeoffset))->second;
            unsigned pinId = edge_id_to_pin[netIdx][nodeoffset];
            const Span& fPt = pinLocations[pinId].first,
                      & sPt = pinLocations[pinId].second;
            if(fPt.first < sPt.first)
            {
                minY  = min(minY, sPt.first);
                maxY  = max(maxY, fPt.second);
                maxYSpans[(*node)->getIndex()].first  = min(maxYSpans[(*node)->getIndex()].first,  fPt.first);
                maxYSpans[(*node)->getIndex()].second = max(maxYSpans[(*node)->getIndex()].second, sPt.second);
            }
            else
            {
                minY  = min(minY, fPt.first);
                maxY  = max(maxY, sPt.second);
                maxYSpans[(*node)->getIndex()].first  = min(maxYSpans[(*node)->getIndex()].first,  sPt.first);
                maxYSpans[(*node)->getIndex()].second = max(maxYSpans[(*node)->getIndex()].second, fPt.second);
            }
        }
        minPossibleYSpan[netIdx] = Span(minY,maxY);
    }

    // Build needUpdate list
    itHGFNodeLocal lc;
    for(n = rbplace.getHGraph().edgesBegin(); n != rbplace.getHGraph().edgesEnd(); ++n)
    {
        unsigned netIdx = (*n)->getIndex();
        Span nSpan = minPossibleYSpan[netIdx];
                   // looking for dependency between nets end cells
                   // this is for checking whether when a cell flipped,
                   // an update to the net is necessary or not
        double netCost = netCosts[netIdx];
        for(lc = (*n)->nodesBegin(); lc != (*n)->nodesEnd(); lc++)
        {
            unsigned cellIdx = (*lc)->getIndex();
            abkfatal(cellIdx < oldCost.size(), "Index too much");
            if(!flipableCells.isBitSet(cellIdx)) continue;
            oldCost[cellIdx] += netCost;
            if(
#ifdef USEFLUTE
               !(_params.useSteiner && (*n)->getDegree() > 3 && (*n)->getDegree() <= steinerCutoff) &&
#endif
               nSpan.first <= maxYSpans[cellIdx].first &&
               maxYSpans[cellIdx].second <= nSpan.second) continue;
            needUpdate[cellIdx].push_back(netIdx);
        }
    }

    // Build gains
    for(c = rbplace.getHGraph().nodesBegin(); c != rbplace.getHGraph().nodesEnd(); ++c)
    {
        unsigned cellIdx = (*c)->getIndex();
        if (!flipableCells.isBitSet(cellIdx)) continue;
        newCost[cellIdx] = oldCost[cellIdx];
        for(unsigned i = 0; i < needUpdate[cellIdx].size(); ++i)
        {
            unsigned netIdx = needUpdate[cellIdx][i];
            newCost[cellIdx] -= netCosts[netIdx];
            const HGFEdge& net = rbplace.getHGraph().getEdgeByIdx(netIdx);
            BBox nBox;
            nBox.xMin = netBBox[netIdx].xMin;
            nBox.xMax = netBBox[netIdx].xMax;
            for(nodeoffset = 0, node = net.nodesBegin(); node != net.nodesEnd(); ++nodeoffset, ++node)
            {
              unsigned pinId = edge_id_to_pin[netIdx][nodeoffset];
              if((*node)->getIndex() == cellIdx) {
                  nBox.yMin = min(nBox.yMin,pinLocations[pinId].second.first);
                  nBox.yMax = max(nBox.yMax,pinLocations[pinId].second.second);
              }
              else {
                  nBox.yMin = min(nBox.yMin,pinLocations[pinId].first.first);
                  nBox.yMax = max(nBox.yMax,pinLocations[pinId].first.second);
              }
            }
            ors[cellIdx].flipVert();
            newCost[cellIdx] += getNetCost(rbplace, netIdx, nBox);
            ors[cellIdx].flipVert();
        }
        if(_params.minimizeCost)
            gain[cellIdx] = oldCost[cellIdx] - newCost[cellIdx];
        else
            gain[cellIdx] = newCost[cellIdx] - oldCost[cellIdx];
    }
    BitBoard visitedCells(numCells);
    HeapGainContainer  heapGain(gain, flipableCells);

//    cout << "gain cutoff is " << _params.cutOffGain << endl;
//    cout << "starting HPWL is " << _initialWL << endl;

                       // Optimization loop starts here
    while(gain[heapGain.getTop()] > _params.cutOffGain)
    {
        visitedCells.clear();
        unsigned cellIdx = heapGain.getTop();
        ors[cellIdx].flipVert();
        _nFlips++;

//        cout << "flip number " << _nFlips
//             << ": best gain is " << gain[heapGain.getTop()]
//             << " for cell " << heapGain.getTop() << endl;

        // problem with swap so use "old-style" one
        double swapBucket = oldCost[cellIdx];
        oldCost[cellIdx] = newCost[cellIdx];
        newCost[cellIdx] = swapBucket;
        gain[cellIdx] = -gain[cellIdx];
        heapGain.heapify(cellIdx);

        if(_params.verb.getForActions() > 1)
            cout << "delta : " << oldCost[cellIdx] - newCost[cellIdx] << endl;
        for(unsigned i = 0; i < needUpdate[cellIdx].size(); ++i)
        {
            unsigned netIdx = needUpdate[cellIdx][i];
            const HGFEdge& net = rbplace.getHGraph().getEdgeByIdx(netIdx);
            BBox nBox;
            for(nodeoffset = 0, node = net.nodesBegin(); node != net.nodesEnd(); ++nodeoffset, ++node) {
                BBox pinOffset = rbplace.getHGraph().nodesOnEdgePinOffsetBBox(nodeoffset,netIdx,ors[(*node)->getIndex()]);
                pinOffset.TranslateBy(rbplace[(*node)->getIndex()]);
                nBox += pinOffset.getBottomLeft();
                nBox += pinOffset.getTopRight();
            }
            double newNetCost = getNetCost(rbplace, netIdx, nBox);
            double oldNetCost = netCosts[netIdx];
            double delta = newNetCost - oldNetCost;
            netBBox[netIdx] = nBox;
            netCosts[netIdx] = newNetCost;
            for(lc = net.nodesBegin(); lc != net.nodesEnd(); lc++)
            {
                unsigned lCellIdx = (*lc)->getIndex();
                if(!flipableCells.isBitSet(lCellIdx)) continue;
                if(lCellIdx == cellIdx) continue;
                oldCost[lCellIdx] += delta;
                if(visitedCells.isBitSet(lCellIdx)) continue;
                visitedCells.setBit(lCellIdx);
            }
        }
        const vector<unsigned>&
            updateCellList = visitedCells.getIndicesOfSetBits();

        for(unsigned i = 0; i < updateCellList.size(); i++)
        {
            unsigned lCellIdx = updateCellList[i];
            ors[lCellIdx].flipVert();
            newCost[lCellIdx] = oldCost[lCellIdx];
            for(unsigned j = 0; j < needUpdate[lCellIdx].size(); j++)
            {
                unsigned netIdx = needUpdate[lCellIdx][j];
                const HGFEdge& net = rbplace.getHGraph().getEdgeByIdx(netIdx);
                BBox nBox;
                for(nodeoffset = 0, node = net.nodesBegin(); node != net.nodesEnd(); ++nodeoffset, ++node) {
                    BBox pinOffset = rbplace.getHGraph().nodesOnEdgePinOffsetBBox(nodeoffset,netIdx,ors[(*node)->getIndex()]);
                    pinOffset.TranslateBy(rbplace[(*node)->getIndex()]);
                    nBox += pinOffset.getBottomLeft();
                    nBox += pinOffset.getTopRight();
                }
                double newNetCost = getNetCost(rbplace, netIdx, nBox);
                double oldNetCost = netCosts[netIdx];
                newCost[lCellIdx] += newNetCost - oldNetCost;
            }
            ors[lCellIdx].flipVert();
            if(_params.minimizeCost)
                gain[lCellIdx] = oldCost[lCellIdx] - newCost[lCellIdx];
            else
                gain[lCellIdx] = newCost[lCellIdx] - oldCost[lCellIdx];
            heapGain.heapify(lCellIdx);
        }
    }
}

void GreedyOrientOpt::optimizeHoriz(const RBPlacement& rbplace, vector<BBox>& netBBox, vector<double> &netCosts,
BitBoard& flipableCells, const vector< vector<unsigned> >& edge_id_to_pin)
{
    unsigned numNets = rbplace.getHGraph().getNumEdges(),
             numCells = rbplace.getHGraph().getNumNodes();
    vector<double> oldCost(numCells, 0.);
    vector<double> newCost(numCells, 0.);
    vector<double> gain(numCells, -DBL_MAX);
    vector<Span> minPossibleXSpan(numNets,Span(DBL_MAX,-DBL_MAX));
    vector<Span> maxXSpans(numCells,Span(DBL_MAX,-DBL_MAX));
    vector< vector<unsigned> > needUpdate(numCells);
    vector<Orient>& ors = _orients;
    vector< pair<Span, Span> > pinLocations(rbplace.getHGraph().getNumPins());

    itHGFNodeGlobal c;
    for(c = rbplace.getHGraph().nodesBegin(); c != rbplace.getHGraph().nodesEnd(); c++)
    {
        if(!rbplace.isCoreCell((*c)->getIndex())) { 
          continue;    // only consider core cells
        }
        if(rbplace.isFixed((*c)->getIndex())) {
          continue; 
        }
        Symmetry symm = rbplace.getHGraph().getNodeSymmetry((*c)->getIndex());
        if(_params.checkSymmetry && !symm.y) {
            continue;        // ignore the cell if it is not Y symmetry
        }
        flipableCells.setBit((*c)->getIndex());
    }
    if(flipableCells.getNumBitsSet() == 0)
    {
        abkwarn(0, "no horiz orientation change is permitted");
        return;
    }

    // Store the actual location of pins for flipped and unflipped
    unsigned edgeoffset, nodeoffset;
    itHGFEdgeLocal e;
    itHGFNodeLocal node;
    unsigned pin_id = 0;
    for(c = rbplace.getHGraph().nodesBegin(); c != rbplace.getHGraph().nodesEnd(); ++c)
    {
        unsigned cellIdx = (*c)->getIndex();
        for(edgeoffset = 0, e = (*c)->edgesBegin(); e != (*c)->edgesEnd(); ++edgeoffset, ++e)
        {
            BBox pinOffset = rbplace.getHGraph().edgesOnNodePinOffsetBBox(edgeoffset,cellIdx,ors[cellIdx]);
            pinOffset.TranslateBy(rbplace[cellIdx]);
            pinLocations[pin_id].first.first  = pinOffset.xMin;
            pinLocations[pin_id].first.second = pinOffset.xMax;
            if(flipableCells.isBitSet(cellIdx))
            {
              ors[cellIdx].flipHoriz();
              pinOffset = rbplace.getHGraph().edgesOnNodePinOffsetBBox(edgeoffset,cellIdx,ors[cellIdx]);
              pinOffset.TranslateBy(rbplace[cellIdx]);
              pinLocations[pin_id].second.first  = pinOffset.xMin;
              pinLocations[pin_id].second.second = pinOffset.xMax;
              ors[cellIdx].flipHoriz();
            }
            else
            {
              pinLocations[pin_id].second = pinLocations[pin_id].first;
            }
            ++pin_id;
        }
    }

    // Store minimum and maximum possible XSpans
    itHGFEdgeGlobal n;
    for(n = rbplace.getHGraph().edgesBegin(); n != rbplace.getHGraph().edgesEnd(); ++n)
    {
        unsigned netIdx = (*n)->getIndex();
        double minX = DBL_MAX, maxX = -DBL_MAX;
        for(nodeoffset = 0, node = (*n)->nodesBegin(); node != (*n)->nodesEnd(); ++nodeoffset, ++node)
        {
            //unsigned pinId = edge_id_to_pin.find(pair<unsigned,unsigned>(netIdx,nodeoffset))->second;
            unsigned pinId = edge_id_to_pin[netIdx][nodeoffset];
            const Span& fPt = pinLocations[pinId].first,
                      & sPt = pinLocations[pinId].second;
            if(fPt.first < sPt.first)
            {
                minX  = min(minX, sPt.first);
                maxX  = max(maxX, fPt.second);
                maxXSpans[(*node)->getIndex()].first  = min(maxXSpans[(*node)->getIndex()].first,  fPt.first);
                maxXSpans[(*node)->getIndex()].second = max(maxXSpans[(*node)->getIndex()].second, sPt.second);
            }
            else
            {
                minX  = min(minX, fPt.first);
                maxX  = max(maxX, sPt.second);
                maxXSpans[(*node)->getIndex()].first  = min(maxXSpans[(*node)->getIndex()].first,  sPt.first);
                maxXSpans[(*node)->getIndex()].second = max(maxXSpans[(*node)->getIndex()].second, fPt.second);
            }
        }
        minPossibleXSpan[netIdx] = Span(minX,maxX);
    }

    // Build needUpdate list
    itHGFNodeLocal lc;
    for(n = rbplace.getHGraph().edgesBegin(); n != rbplace.getHGraph().edgesEnd(); ++n)
    {
        unsigned netIdx = (*n)->getIndex();
        Span nSpan = minPossibleXSpan[netIdx];
                   // looking for dependency between nets end cells
                   // this is for checking whether when a cell flipped,
                   // an update to the net is necessary or not
        double netCost = netCosts[netIdx];
        for(lc = (*n)->nodesBegin(); lc != (*n)->nodesEnd(); lc++)
        {
            unsigned cellIdx = (*lc)->getIndex();
            abkfatal(cellIdx < oldCost.size(), "Index too much");
            if(!flipableCells.isBitSet(cellIdx)) continue;
            oldCost[cellIdx] += netCost;
            if(
#ifdef USEFLUTE
               !(_params.useSteiner && (*n)->getDegree() > 3 && (*n)->getDegree() <= steinerCutoff) &&
#endif
               nSpan.first <= maxXSpans[cellIdx].first &&
               maxXSpans[cellIdx].second <= nSpan.second) continue;
            needUpdate[cellIdx].push_back(netIdx);
        }
    }

    // Build gains
    for(c = rbplace.getHGraph().nodesBegin(); c != rbplace.getHGraph().nodesEnd(); ++c)
    {
        unsigned cellIdx = (*c)->getIndex();
        if (!flipableCells.isBitSet(cellIdx)) continue;
        newCost[cellIdx] = oldCost[cellIdx];
        for(unsigned i = 0; i < needUpdate[cellIdx].size(); ++i)
        {
            unsigned netIdx = needUpdate[cellIdx][i];
            newCost[cellIdx] -= netCosts[netIdx];
            const HGFEdge& net = rbplace.getHGraph().getEdgeByIdx(netIdx);
            BBox nBox;
            nBox.yMin = netBBox[netIdx].yMin;
            nBox.yMax = netBBox[netIdx].yMax;
            for(nodeoffset = 0, node = net.nodesBegin(); node != net.nodesEnd(); ++nodeoffset, ++node)
            {
              unsigned pinId = edge_id_to_pin[netIdx][nodeoffset];
              if((*node)->getIndex() == cellIdx) {
                  nBox.xMin = min(nBox.xMin,pinLocations[pinId].second.first);
                  nBox.xMax = max(nBox.xMax,pinLocations[pinId].second.second);
              }
              else {
                  nBox.xMin = min(nBox.xMin,pinLocations[pinId].first.first);
                  nBox.xMax = max(nBox.xMax,pinLocations[pinId].first.second);
              }
            }
            ors[cellIdx].flipHoriz();
            newCost[cellIdx] += getNetCost(rbplace, netIdx, nBox);
            ors[cellIdx].flipHoriz();
        }
        if(_params.minimizeCost)
            gain[cellIdx] = oldCost[cellIdx] - newCost[cellIdx];
        else
            gain[cellIdx] = newCost[cellIdx] - oldCost[cellIdx];
    }
    BitBoard visitedCells(numCells);
    HeapGainContainer  heapGain(gain, flipableCells);

//    cout << "gain cutoff is " << _params.cutOffGain << endl;
//    cout << "starting HPWL is " << _initialWL << endl;

                       // Optimization loop starts here
    while(gain[heapGain.getTop()] > _params.cutOffGain)
    {
        visitedCells.clear();
        unsigned cellIdx = heapGain.getTop();
        ors[cellIdx].flipHoriz();
        _nFlips++;
  
//        cout << "flip number " << _nFlips
//             << ": best gain is " << gain[heapGain.getTop()]
//             << " for cell " << heapGain.getTop() << endl;

        // problem with swap so use "old-style" one
        double swapBucket = oldCost[cellIdx];
        oldCost[cellIdx] = newCost[cellIdx];
        newCost[cellIdx] = swapBucket;
        gain[cellIdx] = -gain[cellIdx];
        heapGain.heapify(cellIdx);

        if(_params.verb.getForActions() > 1)
            cout << "delta : " << oldCost[cellIdx] - newCost[cellIdx] << endl;
        for(unsigned i = 0; i < needUpdate[cellIdx].size(); ++i)
        {
            unsigned netIdx = needUpdate[cellIdx][i];
            const HGFEdge& net = rbplace.getHGraph().getEdgeByIdx(netIdx);
            BBox nBox;
            for(nodeoffset = 0, node = net.nodesBegin(); node != net.nodesEnd(); ++nodeoffset, ++node) {
                BBox pinOffset = rbplace.getHGraph().nodesOnEdgePinOffsetBBox(nodeoffset,netIdx,ors[(*node)->getIndex()]);
                pinOffset.TranslateBy(rbplace[(*node)->getIndex()]);
                nBox += pinOffset.getBottomLeft();
                nBox += pinOffset.getTopRight();
            }
            double newNetCost = getNetCost(rbplace, netIdx, nBox);
            double oldNetCost = netCosts[netIdx];
            double delta = newNetCost - oldNetCost;
            netBBox[netIdx] = nBox;
            netCosts[netIdx] = newNetCost;
            for(lc = net.nodesBegin(); lc != net.nodesEnd(); lc++)
            {
                unsigned lCellIdx = (*lc)->getIndex();
                if(!flipableCells.isBitSet(lCellIdx)) continue;
                if(lCellIdx == cellIdx) continue;
                oldCost[lCellIdx] += delta;
                if(visitedCells.isBitSet(lCellIdx)) continue;
                visitedCells.setBit(lCellIdx);
            }
        }
        const vector<unsigned>&
            updateCellList = visitedCells.getIndicesOfSetBits();

        for(unsigned i = 0; i < updateCellList.size(); i++)
        {
            unsigned lCellIdx = updateCellList[i];
            ors[lCellIdx].flipHoriz();
            newCost[lCellIdx] = oldCost[lCellIdx];
            for(unsigned j = 0; j < needUpdate[lCellIdx].size(); j++)
            {
                unsigned netIdx = needUpdate[lCellIdx][j];
                const HGFEdge& net = rbplace.getHGraph().getEdgeByIdx(netIdx);
                BBox nBox;
                for(nodeoffset = 0, node = net.nodesBegin(); node != net.nodesEnd(); ++nodeoffset, ++node) {
                    BBox pinOffset = rbplace.getHGraph().nodesOnEdgePinOffsetBBox(nodeoffset,netIdx,ors[(*node)->getIndex()]);
                    pinOffset.TranslateBy(rbplace[(*node)->getIndex()]);
                    nBox += pinOffset.getBottomLeft();
                    nBox += pinOffset.getTopRight();
                }
                double newNetCost = getNetCost(rbplace, netIdx, nBox);
                double oldNetCost = netCosts[netIdx];
                newCost[lCellIdx] += newNetCost - oldNetCost;
            }
            ors[lCellIdx].flipHoriz();
            if(_params.minimizeCost)
                gain[lCellIdx] = oldCost[lCellIdx] - newCost[lCellIdx];
            else
                gain[lCellIdx] = newCost[lCellIdx] - oldCost[lCellIdx];
            heapGain.heapify(lCellIdx);
        }
    }
}
