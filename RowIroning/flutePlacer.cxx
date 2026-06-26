/**************************************************************************
***
*** Copyright (c) 2026 IonQ, Inc.
*** Copyright (c) 1995-2000 Regents of the University of California,
***               Andrew E. Caldwell, Andrew B. Kahng and Igor L. Markov
*** Copyright (c) 2000-2006 Regents of the University of Michigan,
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
#include "ABKCommon/abkcommon.h"
#include "ABKCommon/abklimits.h"
#include "HGraphWDims/hgWDims.h"
#include "RBPlace/RBPlacement.h"
#include "SmallPlacement/smallPlaceProb.h"
#include "Geoms/ptsetpt.h"
#include "rowIroning.h"
#include <map>

using std::vector;

#ifdef USEFLUTE

#include "Flute/flute.h"

double calcLengthOfSize2PtSet(const vector<PtSetPoint> &ptSet)
{
  double netWL  = fabs(ptSet[0].pt.x-ptSet[1].pt.x);
         netWL += fabs(ptSet[0].pt.y-ptSet[1].pt.y);

  return netWL;
}


double calcLengthOfSize3PtSet(const vector<PtSetPoint> &ptSet)
{
    double minX, maxX;

    if(ptSet[0].pt.x < ptSet[1].pt.x)
    {
      minX = ptSet[0].pt.x;
      maxX = ptSet[1].pt.x;
    }
    else
    {
      minX = ptSet[1].pt.x;
      maxX = ptSet[0].pt.x;
    }
    if(ptSet[2].pt.x < minX)
    {
      minX = ptSet[2].pt.x;
    }
    else if(ptSet[2].pt.x > maxX)
    {
      maxX = ptSet[2].pt.x;
    }

    double netWL = (maxX - minX);

    double minY, maxY;

    if(ptSet[0].pt.y < ptSet[1].pt.y)
    {
      minY = ptSet[0].pt.y;
      maxY = ptSet[1].pt.y;
    }
    else
    {
      minY = ptSet[1].pt.y;
      maxY = ptSet[0].pt.y;
    }
    if(ptSet[2].pt.y < minY)
    {
      minY = ptSet[2].pt.y;
    }
    else if(ptSet[2].pt.y > maxY)
    {
      maxY = ptSet[2].pt.y;
    }

    netWL += (maxY - minY);

    return netWL;
}

double calcLengthOfFlutePtSet(const vector<PtSetPoint> &ptSet, double *x, double *y)
{
    for(unsigned p = 0; p < ptSet.size(); ++p)
    {
      x[p] = ptSet[p].pt.x;
      y[p] = ptSet[p].pt.y;
    }

    Tree fltree = flute(ptSet.size(), x, y, 5);
    double netWL = fltree.length;

    free(fltree.branch);

    return netWL;
}

double calcLengthOfLargePtSet(const vector<PtSetPoint> &ptSet)
{
    double xmin = ptSet[0].pt.x;
    double xmax = ptSet.back().pt.x;

    double netWL = xmax - xmin;

    double ymin = ptSet[0].pt.y;
    double ymax = ymin;

    for(unsigned i = 1; i < ptSet.size(); ++i)
    {
      if(ptSet[i].pt.y < ymin)
      {
        ymin = ptSet[i].pt.y;
      }
      else if(ptSet[i].pt.y > ymax)
      {
        ymax = ptSet[i].pt.y;
      }
    }

    netWL += (ymax - ymin);

    return netWL;
}

double calcLengthOfPtSet(const vector<PtSetPoint> &ptSet, double *x, double *y)
{
  double netWL = 0;
  if(ptSet.size() <= 1)
  {
    netWL = 0.;
  }
  else if(ptSet.size() == 2)
  {
    netWL = calcLengthOfSize2PtSet(ptSet);
  }
  else if(ptSet.size() == 3)
  {
    netWL = calcLengthOfSize3PtSet(ptSet);
  }
  else if(ptSet.size() <= 20)
  {
    netWL = calcLengthOfFlutePtSet(ptSet, x, y);
  }
  else
  {
    netWL = calcLengthOfLargePtSet(ptSet);
  }

  return netWL;
}

void addToPointSet(vector<PtSetPoint> &ptSet, const Point &pt)
{
  PtSetPoint ptSetPt(pt, 1);

  std::pair<vector<PtSetPoint>::iterator, vector<PtSetPoint>::iterator> iters =
       std::equal_range(ptSet.begin(), ptSet.end(), ptSetPt);

  if(iters.first == iters.second)
  {
    ptSet.insert(iters.first, ptSetPt);
  }
  else
  {
    abkfatal((iters.first+1) == iters.second, "duplicates found in point set");
    abkfatal(iters.first->count >= 1, "bad point set point");
    ++iters.first->count;
  }
}

void removeFromPointSet(vector<PtSetPoint> &ptSet, const Point &pt)
{
  PtSetPoint ptSetPt(pt, 1);

  std::pair<vector<PtSetPoint>::iterator, vector<PtSetPoint>::iterator> iters =
       std::equal_range(ptSet.begin(), ptSet.end(), ptSetPt);

  abkfatal(iters.first != iters.second, "pt to remove not found in point set");
  abkfatal((iters.first+1) == iters.second, "duplicates found in point set");
  abkfatal(iters.first->count >= 1, "bad point set point");

  if(iters.first->count > 1)
  {
    --iters.first->count;
  }
  else
  {
    ptSet.erase(iters.first);
  }
}

void placeWithFlute(const RBPlacement &rbplace,
                    const vector< vector<unsigned> > &movables,
                    const SmallPlacementRow &row,
                    const vector<double> &whiteSpaceWidths,
                    BitBoard &nets, Placement &soln, Placement &currSoln,
                    vector< vector<PtSetPoint> > &ptSets)
{
   const HGraphWDimensions &hgWDims = rbplace.getHGraph();

   nets.clear();

   for(unsigned i = 0; i < movables.size(); ++i)
   {
     if(movables[i][0] == UINT_MAX) continue;

     for(unsigned j = 0; j < movables[i].size(); ++j)
     {
       const HGFNode &node = hgWDims.getNodeByIdx(movables[i][j]);

       for(itHGFEdgeLocal e = node.edgesBegin(); e != node.edgesEnd(); ++e)
       {
         nets.setBit((*e)->getIndex());
       }
     }
   }

   const vector<unsigned>& affectedNets = nets.getIndicesOfSetBits();
   const unsigned numNets = affectedNets.size();

   vector<unsigned> ordering(movables.size());
   vector<double> cellWidths(movables.size(),0.);

   unsigned whiteSpaceCell = 0;
   for(unsigned i = 0; i < movables.size(); ++i)
   {
     ordering[i] = i;
     if(movables[i][0] == UINT_MAX)
     {
       cellWidths[i] = whiteSpaceWidths[whiteSpaceCell++];
     }
     else
     {
       for(unsigned j = 0; j < movables[i].size(); ++j)
       {
         cellWidths[i] += hgWDims.getNodeWidth(movables[i][j]);
       }
     }
   }

   double x[20], y[20];
   double bestWL = std::numeric_limits<double>::max();

   do
   {
     bool skip_this_perm = false;
     for(unsigned i = 1; i < movables.size(); ++i)
     {
        if(movables[ordering[i-1]][0] == UINT_MAX  &&
           movables[ordering[i]][0] == UINT_MAX &&
           ordering[i-1] < ordering[i])
        {
           // avoid orientations that will end with the exact same score
           // i.e. situations where just fake cells have been swapped
           skip_this_perm = true;
           break;
        }
     }

     if(skip_this_perm)
     {
       continue;
     }

     // calculate the x coords of the cells
     Point pos = row.getBottomLeft();
     for(unsigned i = 0; i < movables.size(); ++i)
     {
        if(movables[ordering[i]][0] != UINT_MAX)
        {
          if(movables[ordering[i]].size() == 1)
          {
            unsigned cellId = movables[ordering[i]][0];

            // update the point sets
            if(currSoln[cellId] != pos)
            {
              const HGFNode &cell = hgWDims.getNodeByIdx(cellId);
              unsigned offset = 0;
              for(itHGFEdgeLocal e = cell.edgesBegin(); e != cell.edgesEnd(); ++e, ++offset)
              {
                unsigned netId = (*e)->getIndex();
                BBox pinOffset = hgWDims.edgesOnNodePinOffsetBBox(offset, cellId, rbplace.getOrient(cellId));

                removeFromPointSet(ptSets[netId], currSoln[cellId] + pinOffset.getBottomLeft());
                addToPointSet(ptSets[netId], pos + pinOffset.getBottomLeft());

                if(!equalDouble(pinOffset.xMin,pinOffset.xMax) ||
                   !equalDouble(pinOffset.yMin,pinOffset.yMax))
                {
                  removeFromPointSet(ptSets[netId], currSoln[cellId] + pinOffset.getTopRight());
                  addToPointSet(ptSets[netId], pos + pinOffset.getTopRight());
                }
              }

              currSoln[cellId] = pos;
            }
            pos.x += cellWidths[ordering[i]];
          }
          else
          {
            for(unsigned j = 0; j < movables[ordering[i]].size(); ++j)
            {
              unsigned cellId = movables[ordering[i]][j];

              // update the point sets
              if(currSoln[cellId] != pos)
              {
                const HGFNode &cell = hgWDims.getNodeByIdx(cellId);
                unsigned offset = 0;
                for(itHGFEdgeLocal e = cell.edgesBegin(); e != cell.edgesEnd(); ++e, ++offset)
                {
                  unsigned netId = (*e)->getIndex();
                  BBox pinOffset = hgWDims.edgesOnNodePinOffsetBBox(offset, cellId, rbplace.getOrient(cellId));

                  removeFromPointSet(ptSets[netId], currSoln[cellId] + pinOffset.getBottomLeft());
                  addToPointSet(ptSets[netId], pos + pinOffset.getBottomLeft());

                  if(!equalDouble(pinOffset.xMin,pinOffset.xMax) ||
                     !equalDouble(pinOffset.yMin,pinOffset.yMax))
                  {
                    removeFromPointSet(ptSets[netId], currSoln[cellId] + pinOffset.getTopRight());
                    addToPointSet(ptSets[netId], pos + pinOffset.getTopRight());
                  }
                }

                currSoln[cellId] = pos;
              }
              pos.x += hgWDims.getNodeWidth(cellId);
            }
          }
        }
        else
        {
          pos.x += cellWidths[ordering[i]];
        }
     }

     // calculate the cost of this ordering
     double curWL = 0.;
     for(unsigned i = 0; i < numNets; ++i)
     {
       curWL += calcLengthOfPtSet(ptSets[affectedNets[i]], x, y);
     }

     if(curWL < bestWL)
     {
        bestWL = curWL;
        for(unsigned i = 0; i < movables.size(); ++i)
        {
          if(movables[i][0] == UINT_MAX) continue;
          soln[i] = currSoln[movables[i][0]];
        }
     }
   }
   while(next_permutation(ordering.begin(), ordering.end()));
}

class DPTableEntry
{
public:
   Placement place;
   std::map<unsigned,double> netLens;
   Point rightEdge;
   double cost;

   DPTableEntry(unsigned numCells, Point edge):
                place(numCells), rightEdge(edge), cost(0.) 
   {
     for(unsigned i = 0; i < numCells; ++i)
     {
       place[i].x = place[i].y = DBL_MAX;
     }
   }
};

void placeWithFluteDP(const RBPlacement &rbplace,
                      const vector< vector<unsigned> > &movables,
                      const SmallPlacementRow &row,
                      const vector<double> &whiteSpaceWidths,
                      BitBoard &nets, Placement &soln, Placement &tempSoln,
                      vector< vector<PtSetPoint> > &ptSets)
{
   const HGraphWDimensions &hgWDims = rbplace.getHGraph();
   std::map<unsigned,double> tempNetLensA, tempNetLensB;
   vector<double> cellWidths(movables.size(),0.);

   unsigned whiteSpaceCell = 0;
   for(unsigned i = 0; i < movables.size(); ++i)
   {
     if(movables[i][0] == UINT_MAX)
     {
       cellWidths[i] = whiteSpaceWidths[whiteSpaceCell++];
     }
     else
     {
       for(unsigned j = 0; j < movables[i].size(); ++j)
       {
         cellWidths[i] += hgWDims.getNodeWidth(movables[i][j]);
       }
     }
   }

   // calculate the initial cost, we don't want to be worse
   double x[20], y[20];
   double startWL = 0.;
   nets.clear();
   for(unsigned i = 0; i < movables.size(); ++i)
   {
     if(movables[i][0] == UINT_MAX) continue;
     for(unsigned j = 0; j < movables[i].size(); ++j)
     {
       const HGFNode &node = hgWDims.getNodeByIdx(movables[i][j]);
       for(itHGFEdgeLocal e = node.edgesBegin(); e != node.edgesEnd(); ++e)
       {
         unsigned netId = (*e)->getIndex();
         if(nets.isBitSet(netId)) continue;
         nets.setBit(netId);
         startWL += calcLengthOfPtSet(ptSets[netId], x, y);
       }
     }
   }

   unsigned numAcells = movables.size()/2;
   unsigned numBcells = movables.size()-numAcells;
   vector< vector< DPTableEntry > > dpTable(numAcells+1,
           vector< DPTableEntry >(numBcells+1, DPTableEntry(movables.size(),row.getBottomLeft())));

   Point notPlaced(DBL_MAX,DBL_MAX);

   // fill in the entries of the table
   for(unsigned i = 0; i <= numAcells; ++i)
   {
     for(unsigned j = 0; j <= numBcells; ++j)
     {
       double costA = DBL_MAX;
       if(i > 0)
       {
         costA = dpTable[i-1][j].cost;

         //copy net lengths into temp vector
         tempNetLensA = dpTable[i-1][j].netLens;

         if(movables[i-1][0] != UINT_MAX)
         {
           //copy locations of dpTable[i-1][j] into the temp placement
           for(unsigned k = 0; k < movables.size(); ++k)
           {
             if(movables[k][0] != UINT_MAX)
             {
               Point pos = dpTable[i-1][j].place[k];
               for(unsigned l = 0; l < movables[k].size(); ++l)
               {
                  const unsigned cellId = movables[k][l];
                  // update the point sets
                  if(tempSoln[cellId] != pos)
                  {
                    const HGFNode &cell = hgWDims.getNodeByIdx(cellId);
                    unsigned offset = 0;
                    for(itHGFEdgeLocal e = cell.edgesBegin(); e != cell.edgesEnd(); ++e, ++offset)
                    {
                      unsigned netId = (*e)->getIndex();
                      BBox pinOffset = hgWDims.edgesOnNodePinOffsetBBox(offset, cellId, rbplace.getOrient(cellId));

                      if(tempSoln[cellId] != notPlaced)
                         removeFromPointSet(ptSets[netId], tempSoln[cellId] + pinOffset.getBottomLeft());
                      if(pos != notPlaced)
                         addToPointSet(ptSets[netId], pos + pinOffset.getBottomLeft());

                      if(!equalDouble(pinOffset.xMin,pinOffset.xMax) ||
                         !equalDouble(pinOffset.yMin,pinOffset.yMax))
                      {
                        if(tempSoln[cellId] != notPlaced)
                          removeFromPointSet(ptSets[netId], tempSoln[cellId] + pinOffset.getTopRight());
                        if(pos != notPlaced)
                          addToPointSet(ptSets[netId], pos + pinOffset.getTopRight());
                      }
                   }

                   tempSoln[cellId] = pos;
                 }
                 pos.x += hgWDims.getNodeWidth(cellId);
               }
             }
           }

           //place the new cell(s)
           Point pos = dpTable[i-1][j].rightEdge;
           for(unsigned l = 0; l < movables[i-1].size(); ++l)
           {
             const unsigned cellId = movables[i-1][l];
             // update the point sets
             if(tempSoln[cellId] != pos)
             {
               const HGFNode &cell = hgWDims.getNodeByIdx(cellId);
               unsigned offset = 0;
               for(itHGFEdgeLocal e = cell.edgesBegin(); e != cell.edgesEnd(); ++e, ++offset)
               {
                 unsigned netId = (*e)->getIndex();
                 BBox pinOffset = hgWDims.edgesOnNodePinOffsetBBox(offset, cellId, rbplace.getOrient(cellId));

                 if(tempSoln[cellId] != notPlaced)
                   removeFromPointSet(ptSets[netId], tempSoln[cellId] + pinOffset.getBottomLeft());
                 addToPointSet(ptSets[netId], pos + pinOffset.getBottomLeft());

                 if(!equalDouble(pinOffset.xMin,pinOffset.xMax) ||
                    !equalDouble(pinOffset.yMin,pinOffset.yMax))
                 {
                   if(tempSoln[cellId] != notPlaced)
                     removeFromPointSet(ptSets[netId], tempSoln[cellId] + pinOffset.getTopRight());
                   addToPointSet(ptSets[netId], pos + pinOffset.getTopRight());
                 }
              }

              tempSoln[cellId] = pos;
             }
             pos.x += hgWDims.getNodeWidth(cellId);
           }

           //figure out the new cost
           nets.clear();
           for(unsigned l = 0; l < movables[i-1].size(); ++l)
           {
             const HGFNode &node = hgWDims.getNodeByIdx(movables[i-1][l]);
             for(itHGFEdgeLocal e = node.edgesBegin(); e != node.edgesEnd(); ++e)
             {
               unsigned netId = (*e)->getIndex();
               if(nets.isBitSet(netId)) continue;
               nets.setBit(netId);
               std::map<unsigned,double>::const_iterator nl = tempNetLensA.find(netId);
               if(nl != tempNetLensA.end()) { costA -= nl->second; }
               double netLen = calcLengthOfPtSet(ptSets[netId], x, y);
               tempNetLensA[netId] = netLen;
               costA += netLen;
             }
           }
         }
       }

       double costB = DBL_MAX;
       if(j > 0)
       {
         costB = dpTable[i][j-1].cost;

         //copy net lengths into temp vector
         tempNetLensB = dpTable[i][j-1].netLens;

         if(movables[numAcells+j-1][0] != UINT_MAX)
         {
           //copy locations of dpTable[i][j-1] into the temp placement
           for(unsigned k = 0; k < movables.size(); ++k)
           {
             if(movables[k][0] != UINT_MAX)
             {
               Point pos = dpTable[i][j-1].place[k];
               for(unsigned l = 0; l < movables[k].size(); ++l)
               {
                  const unsigned cellId = movables[k][l];
                  // update the point sets
                  if(tempSoln[cellId] != pos)
                  {
                    const HGFNode &cell = hgWDims.getNodeByIdx(cellId);
                    unsigned offset = 0;
                    for(itHGFEdgeLocal e = cell.edgesBegin(); e != cell.edgesEnd(); ++e, ++offset)
                    {
                      unsigned netId = (*e)->getIndex();
                      BBox pinOffset = hgWDims.edgesOnNodePinOffsetBBox(offset, cellId, rbplace.getOrient(cellId));

                      if(tempSoln[cellId] != notPlaced)
                        removeFromPointSet(ptSets[netId], tempSoln[cellId] + pinOffset.getBottomLeft());
                      if(pos != notPlaced)
                      addToPointSet(ptSets[netId], pos + pinOffset.getBottomLeft());

                      if(!equalDouble(pinOffset.xMin,pinOffset.xMax) ||
                         !equalDouble(pinOffset.yMin,pinOffset.yMax))
                      {
                        if(tempSoln[cellId] != notPlaced)
                          removeFromPointSet(ptSets[netId], tempSoln[cellId] + pinOffset.getTopRight());
                        if(pos != notPlaced)
                          addToPointSet(ptSets[netId], pos + pinOffset.getTopRight());
                      }
                   }

                   tempSoln[cellId] = pos;
                 }
                 pos.x += hgWDims.getNodeWidth(cellId);
               }
             }
           }

           //place the new cell(s)
           Point pos = dpTable[i][j-1].rightEdge;
           for(unsigned l = 0; l < movables[numAcells+j-1].size(); ++l)
           {
             const unsigned cellId = movables[numAcells+j-1][l];
             // update the point sets
             if(tempSoln[cellId] != pos)
             {
               const HGFNode &cell = hgWDims.getNodeByIdx(cellId);
               unsigned offset = 0;
               for(itHGFEdgeLocal e = cell.edgesBegin(); e != cell.edgesEnd(); ++e, ++offset)
               {
                 unsigned netId = (*e)->getIndex();
                 BBox pinOffset = hgWDims.edgesOnNodePinOffsetBBox(offset, cellId, rbplace.getOrient(cellId));

                 if(tempSoln[cellId] != notPlaced)
                   removeFromPointSet(ptSets[netId], tempSoln[cellId] + pinOffset.getBottomLeft());
                 addToPointSet(ptSets[netId], pos + pinOffset.getBottomLeft());

                 if(!equalDouble(pinOffset.xMin,pinOffset.xMax) ||
                    !equalDouble(pinOffset.yMin,pinOffset.yMax))
                 {
                   if(tempSoln[cellId] != notPlaced)
                     removeFromPointSet(ptSets[netId], tempSoln[cellId] + pinOffset.getTopRight());
                   addToPointSet(ptSets[netId], pos + pinOffset.getTopRight());
                 }
              }

              tempSoln[cellId] = pos;
             }
             pos.x += hgWDims.getNodeWidth(cellId);
           }

           //figure out the new cost
           nets.clear();
           for(unsigned l = 0; l < movables[numAcells+j-1].size(); ++l)
           {
             const HGFNode &node = hgWDims.getNodeByIdx(movables[numAcells+j-1][l]);
             for(itHGFEdgeLocal e = node.edgesBegin(); e != node.edgesEnd(); ++e)
             {
               unsigned netId = (*e)->getIndex();
               if(nets.isBitSet(netId)) continue;
               nets.setBit(netId);
               std::map<unsigned,double>::const_iterator nl = tempNetLensB.find(netId);
               if(nl != tempNetLensB.end()) { costB -= nl->second; }
               double netLen = calcLengthOfPtSet(ptSets[netId], x, y);
               tempNetLensB[netId] = netLen;
               costB += netLen;
             }
           }
         }
       }

       if(costA < costB)
       {
         dpTable[i][j].place = dpTable[i-1][j].place;
         dpTable[i][j].place[i-1] = dpTable[i-1][j].rightEdge;
         dpTable[i][j].netLens = tempNetLensA;
         dpTable[i][j].rightEdge = dpTable[i-1][j].rightEdge;
         dpTable[i][j].rightEdge.x += cellWidths[i-1];
         dpTable[i][j].cost = costA;
       }
       else if(costB < DBL_MAX)
       {
         dpTable[i][j].place = dpTable[i][j-1].place;
         dpTable[i][j].place[numAcells+j-1] = dpTable[i][j-1].rightEdge;
         dpTable[i][j].netLens = tempNetLensB;
         dpTable[i][j].rightEdge = dpTable[i][j-1].rightEdge;
         dpTable[i][j].rightEdge.x += cellWidths[numAcells+j-1];
         dpTable[i][j].cost = costB;
       }
     }
   }

   // examine table entry dpTable[numAcells][numBcells]
   // if the cost is better than what we started with, use the soln
   if(dpTable[numAcells][numBcells].cost <= startWL)
   {
     for(unsigned i = 0; i < movables.size(); ++i)
     {
       soln[i] = dpTable[numAcells][numBcells].place[i];
     }
   }
   else
   {
     for(unsigned i = 0; i < movables.size(); ++i)
     {
       if(movables[i][0] != UINT_MAX)
       {
         soln[i] = rbplace.getPlacement()[movables[i][0]];
       }
     }
   }
}

#endif
