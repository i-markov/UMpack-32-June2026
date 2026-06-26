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


#include "rowIroning.h"
#include <algorithm>
#include <iostream>
#include "HGraphWDims/hgWDims.h"
#include "RBPlace/RBPlacement.h"
#include "RBPlace/pinLocCalc.h"
#include "SmallPlacement/smallPlaceProb.h"
#include "SmallPlacers/sRowBBPlacer.h"
#include "SmallPlacers/sRowSmPlacer.h"
#include "PlaceEvals/pinPlEval.h"
#include "ABKCommon/abkcommon.h"
#ifdef USEFLUTE
#include "Flute/flute.h"
#endif

using std::cout;
using std::endl;
using std::vector;

void ironRows(RBPlacement& rbplace, SmallPlacementBitBoardContainer &spBitBoards, MaxMem *maxMem)
{
  RowIroningParameters RIParams;
  ironRows(rbplace, RIParams, spBitBoards, maxMem);
}

void ironRows(RBPlacement& rbplace, const RowIroningParameters& params,
              SmallPlacementBitBoardContainer &spBitBoards, MaxMem *maxMem)
{
    Timer totalTime;

    // check that populated ?
    // check params

    abkfatal(params.windowSize>1 && params.windowSize<20,
            " Row ironing window size must be >1 and <20");
    abkfatal(                    params.overlap<params.windowSize,
            " Row ironing overlap must be        smaller than window size\n");
    abkfatal(params.numPasses<20, "Row ironing cannot have 20+ passes");

    if(params.useDP) {
        params.smplParams.algo = SmallPlParams::DynamicProgram;
    }

#ifdef USEFLUTE
    if(params.useSteiner)
    {
      readLUT();
    }
#endif

    double initHPWL = rbplace.evalHPWL(true);
    for(unsigned pass=0; pass<params.numPasses; pass++)
    {
        Timer passTime;
        bool twoDim = 0;
        if(params.twoDim)
            twoDim = 1;

        bool cluster = params.cluster;

        if(params.mixed && !params.twoDim)
        {
            if(pass%2 == 0)
                twoDim = 1;
            else
                twoDim = 0;

            //alternate vertical directions in groups of two
            bool verticalDir=false;
            if( (pass/2)%2 == 0)
                verticalDir=false;
            else
                verticalDir=true;

            //alternate left to right in groups of four
            //starting with that specified in params.LRfirst
            bool leftToRight=false;
            if( (pass/4)%2 == 0 )
                leftToRight=true;
            else
                leftToRight=false;
            if(!params.LRfirst)
                leftToRight=!leftToRight;

            ironRows2D(rbplace, params, spBitBoards, maxMem, verticalDir, twoDim, cluster, leftToRight);
        }
        else
        {
            bool verticalDir=false;
            if(pass%2 == 0)
                verticalDir = false;
            else
                verticalDir = true;

            //alternate left to right in groups of two
            //starting with that specified in params.LRfirst
            bool leftToRight=false;
            if( (pass/2)%2 == 0 )
                leftToRight=true;
            else
                leftToRight=false;
            if(!params.LRfirst)
                leftToRight=!leftToRight;

            ironRows2D(rbplace, params, spBitBoards, maxMem, verticalDir, twoDim, cluster, leftToRight);

        }
        passTime.stop();
        if (params.verb.getForMajStats() || params.verb.getForSysRes() )
        {
            cout << " Pass " << pass << " took : " << passTime << endl;
            cout << "  Pin-to-Pin HPWL after pass: "
                <<rbplace.evalHPWL(true)<<endl;
        }
    }
    totalTime.stop();
    double finalHPWL = rbplace.evalHPWL(true);

    if (params.verb.getForMajStats() || params.verb.getForSysRes() )
    {
        cout << " Row ironing took : " << totalTime << endl;
        double impr = (initHPWL == 0) ? 0 : (initHPWL-finalHPWL)/initHPWL;
        if(impr < 0) cout << "  % increase in HPWL is " << floor(-impr*100000+0.5)/1000<<endl;
        else cout << "  % decrease in HPWL is " << floor(impr*100000+0.5)/1000<<endl;
    }

    //rbplace.snapCellsInSites();
    //rbplace.remOverlaps();
    rbplace.resetPlacement();
}

void ironRows2D(RBPlacement& rbplace, const RowIroningParameters& params,
                SmallPlacementBitBoardContainer &spBitBoards, MaxMem *maxMem, bool verticalDir,
                bool twoDim, bool cluster, bool leftToRight)
{
    if(leftToRight)
		ironRows2DLR(rbplace, params, spBitBoards, maxMem, verticalDir, twoDim, cluster);
	else
		ironRows2DRL(rbplace, params, spBitBoards, maxMem, verticalDir, twoDim, cluster);
}


//to iron in two dimension left to right
void ironRows2DLR(RBPlacement& rbplace, const RowIroningParameters& params,
                  SmallPlacementBitBoardContainer &spBitBoards, MaxMem *maxMem,
                  bool verticalDir, bool twoDim, bool cluster)
{
    const HGraphWDimensions& hgWDims = rbplace.getHGraph();
#ifdef USEFLUTE
    BitBoard forFlute(hgWDims.getNumEdges());
    Placement tempSolnForFlute(rbplace.getPlacement());

    vector< vector<PtSetPoint> > ptSets(hgWDims.getNumEdges());

    // populate the point sets
    for(unsigned i = 0; i < ptSets.size(); ++i)
    {
      const HGFEdge &net = hgWDims.getEdgeByIdx(i);
      ptSets[i].reserve(2*net.getDegree());

      unsigned offset = 0;
      for(itHGFNodeLocal c = net.nodesBegin(); c != net.nodesEnd(); ++c, ++offset)
      {
        unsigned cellId = (*c)->getIndex();
        BBox pinOffset = hgWDims.nodesOnEdgePinOffsetBBox(offset, i, rbplace.getOrient(cellId));

        addToPointSet(ptSets[i], tempSolnForFlute[cellId] + pinOffset.getBottomLeft());

        if(!equalDouble(pinOffset.xMin,pinOffset.xMax) ||
           !equalDouble(pinOffset.yMin,pinOffset.yMax))
        {
          addToPointSet(ptSets[i], tempSolnForFlute[cellId] + pinOffset.getTopRight());
        }
      }
    }
#endif
    bit_vector placed = bit_vector(rbplace.getNumCells(), 1);
    RandomDouble rndDBL(-1, 1);
//    char fileName1[15];
//    char fileName2[15];
//    char fileName3[15];
//    char fileName4[15];
//    unsigned instance = 0;
    unsigned largeNetProblems = 0;

    SmallPlParams::SmPlAlgo algo_save;
    if(twoDim) {
        algo_save = params.smplParams.algo;
        params.smplParams.algo = SmallPlParams::DynamicProgram;
    }  

    //verticalDir=0 means bottom up pass else
    //top down pass
    itRBPCoreRow rIt1, rIt2;
    itRBPSubRow sIt1, sIt2;

    if(verticalDir == 0)
        rIt1=rbplace.coreRowsBegin();
    else
        rIt1=rbplace.coreRowsEnd()-1;

    while(true)                 //exit only by a break
    {
        unsigned windowSize1=params.windowSize;
        unsigned windowSize2=0;
        rIt2 = rIt1;
        if(verticalDir == 0 && twoDim == 1)
        {
            if(rIt1 != rbplace.coreRowsEnd()-1)   //if not last core row
            {
                ++rIt2;
                windowSize1 = unsigned(ceil(params.windowSize/2.0));
                windowSize2 = params.windowSize - windowSize1;
            }
        }
        else
        {
            if(rIt1 != rbplace.coreRowsBegin() && twoDim == 1) 
                //if not first Row
            {
                --rIt2;
                windowSize1 = unsigned(ceil(params.windowSize/2.0));
                windowSize2 = params.windowSize - windowSize1;
            }
        }

        sIt1=rIt1->subRowsBegin();
        sIt2=rIt2->subRowsBegin();

        while(sIt1 != rIt1->subRowsEnd() && sIt2 != rIt2->subRowsEnd())
        {
            const RBPSubRow& subrow1=*sIt1;
            const RBPSubRow& subrow2=*sIt2;

            sIt1++;
            sIt2++;

            SmallPlacementRow smallRow1;
            SmallPlacementRow smallRow2;

            smallRow1.siteInterval=subrow1.getSpacing();
            smallRow1.xMax =smallRow1.xMin=subrow1.getXStart();
            smallRow1.yMax =smallRow1.yMin=subrow1.getYCoord();
            smallRow1.yMax+=subrow1.getHeight();
            unsigned firstCell1=0;
            unsigned totalCells1=subrow1.getNumCells();
            unsigned currCellLoc1=0;
            unsigned nodeCtr1=0;

            smallRow2.siteInterval=subrow2.getSpacing();
            smallRow2.xMax =smallRow2.xMin=subrow2.getXStart();
            smallRow2.yMax =smallRow2.yMin=subrow2.getYCoord();
            smallRow2.yMax+=subrow2.getHeight();
            unsigned firstCell2=0;
            unsigned totalCells2=subrow2.getNumCells();
            unsigned currCellLoc2=0;
            unsigned nodeCtr2=0;

            int numSitesWS=0;
            int sitesPerWS=0;

            if(totalCells1 < 1 || totalCells2 < 1) continue;

            //now start forming instances
            while (firstCell1 < totalCells1)
            {   
                if(totalCells1 < 1 || totalCells2 < 1) break;

//                ++instance;
//                sprintf(fileName1,"smPbIn-%d",instance);
//                sprintf(fileName2,"smPbOut-%d",instance);
//                sprintf(fileName3,"rbPbIn-%d",instance);
//                sprintf(fileName4,"rbPbOut-%d",instance);

                vector<SmallPlacementRow> smallPbRows;
                //form the 1'st row first
                //get the cells + whitespace cells here
                vector< vector<unsigned> > movables;
                vector<double> whiteSpaceWidths;
                vector<Point> whiteSpaceLocs;
                unsigned totalMovables1=0;
                unsigned totalMovables2=0;
                unsigned row1Size=0;

                for(nodeCtr2 = 0; nodeCtr2 < totalCells2-1; ++nodeCtr2)
                {
                    if((rbplace[*(subrow2.cellIdsBegin()+nodeCtr2)].x <=
                        rbplace[*(subrow1.cellIdsBegin()+firstCell1)].x) &&
                       (rbplace[*(subrow2.cellIdsBegin()+nodeCtr2+1)].x >
                        rbplace[*(subrow1.cellIdsBegin()+firstCell1)].x))
                    {
                        break;
                    }
                }
                firstCell2 = nodeCtr2;

                smallRow1.xMax = smallRow1.xMin = rbplace[*(subrow1.cellIdsBegin()+firstCell1)].x;

                for(nodeCtr1 = firstCell1; nodeCtr1 < totalCells1;)
                {
                    vector <unsigned> cellIds;
                    //if not last cell
                    if(nodeCtr1 < totalCells1-1)
                    {
                        unsigned currentCell = *(subrow1.cellIdsBegin()+nodeCtr1);
                        unsigned nextCell = *(subrow1.cellIdsBegin()+nodeCtr1+1);
                        double totalWhiteSpace = rbplace[nextCell].x - rbplace[currentCell].x - hgWDims.getNodeWidth(currentCell);

                        //if whitespace found
                        if(totalWhiteSpace > 0)
                        {
                            cellIds.push_back(currentCell);
                            movables.push_back(cellIds);
                            cellIds.pop_back();
                            smallRow1.xMax = rbplace[currentCell].x + hgWDims.getNodeWidth(currentCell);
                            ++nodeCtr1;
                            ++totalMovables1;
                            if(totalMovables1 == windowSize1)
                              break;

                            Point currWSLoc;
                            currWSLoc.x = smallRow1.xMax;
                            currWSLoc.y = smallRow1.yMin;
                            bool LR = true;
                            bool movablesFull = divideWhitespace(subrow1, windowSize1,
                                 LR, smallRow1, totalWhiteSpace, currWSLoc, totalMovables1, movables,
                                 whiteSpaceWidths, whiteSpaceLocs, params.wsChunkCap);

                            if(movablesFull)
                              break;
                        }
                        else //this cell has no WS in back of it
                        {
                            double rd = rndDBL;

                            cellIds.push_back(currentCell);
                            smallRow1.xMax = rbplace[currentCell].x + hgWDims.getNodeWidth(currentCell);
                            ++nodeCtr1;

                            //cluster the next cell based on random value
                            if(rd > 0 && cluster == 1)   //cluster
                            {
                                cellIds.push_back(nextCell);
                                smallRow1.xMax = rbplace[nextCell].x + hgWDims.getNodeWidth(nextCell);
                                ++nodeCtr1;
                            }
                            movables.push_back(cellIds);
                            cellIds.clear();
                            ++totalMovables1;
                            if(totalMovables1 == windowSize1)
                                break;

                            //if white space present after clustered cell
                            if(nodeCtr1 < totalCells1 && rd > 0 && cluster == 1)
                            {
                                unsigned thirdCell = *(subrow1.cellIdsBegin()+nodeCtr1);
                                totalWhiteSpace = rbplace[thirdCell].x - rbplace[nextCell].x - hgWDims.getNodeWidth(nextCell);
                                Point currWSLoc;
                                if(totalWhiteSpace >= 2*smallRow1.siteInterval)
                                {
                                    numSitesWS = int(totalWhiteSpace/smallRow1.siteInterval);
                                    cellIds.push_back(UINT_MAX);
                                    movables.push_back(cellIds);
                                    cellIds.pop_back();
                                    sitesPerWS = int(floor(numSitesWS/2+0.5));
                                    numSitesWS -= sitesPerWS;
                                    whiteSpaceWidths.push_back(sitesPerWS*smallRow1.siteInterval);
                                    currWSLoc.x=smallRow1.xMax;
                                    currWSLoc.y=smallRow1.yMin;
                                    whiteSpaceLocs.push_back(currWSLoc);
                                    currWSLoc.x += sitesPerWS*smallRow1.siteInterval;
                                    smallRow1.xMax += sitesPerWS*smallRow1.siteInterval;
                                    totalWhiteSpace -= sitesPerWS*smallRow1.siteInterval;
                                    ++totalMovables1;
                                    if(totalMovables1 == windowSize1)
                                        break;

                                    cellIds.push_back(UINT_MAX);
                                    movables.push_back(cellIds);
                                    cellIds.pop_back();
                                    sitesPerWS = numSitesWS;
                                    whiteSpaceWidths.push_back(sitesPerWS*smallRow1.siteInterval);
                                    whiteSpaceLocs.push_back(currWSLoc);
                                    currWSLoc.x += sitesPerWS*smallRow1.siteInterval ;
                                    smallRow1.xMax += sitesPerWS*smallRow1.siteInterval;
                                    totalWhiteSpace -= sitesPerWS*smallRow1.siteInterval;
                                    ++totalMovables1;
                                    if(totalMovables1 == windowSize1)
                                        break;
                                }
                                else if(totalWhiteSpace > 0)
                                {
                                    cellIds.push_back(UINT_MAX);
                                    movables.push_back(cellIds);
                                    cellIds.pop_back();
                                    whiteSpaceWidths.push_back(totalWhiteSpace);
                                    currWSLoc.x=smallRow1.xMax;
                                    currWSLoc.y=smallRow1.yMin;
                                    whiteSpaceLocs.push_back(currWSLoc);
                                    currWSLoc.x += totalWhiteSpace;
                                    smallRow1.xMax += totalWhiteSpace;
                                    ++totalMovables1;
                                    if(totalMovables1 == windowSize1)
                                        break;
                                }
                            }
                        }
                    }
                    else
                    {
                        cellIds.push_back(*(subrow1.cellIdsBegin()+nodeCtr1));
                        movables.push_back(cellIds);
                        cellIds.pop_back();
                        smallRow1.xMax=
                            rbplace[*(subrow1.cellIdsBegin()+nodeCtr1)].x+
                            hgWDims.getNodeWidth(*(subrow1.cellIdsBegin()+nodeCtr1));
                        ++nodeCtr1;
                        ++totalMovables1;
                        if(totalMovables1 == windowSize1)
                            break;		  

                        //if whitespace present after last cell
                        double currWSWidth=subrow1.getXEnd() - smallRow1.xMax;
                        Point currWSLoc;
                        if(currWSWidth >= 2*smallRow1.siteInterval)
                        {
                            numSitesWS = int(currWSWidth/smallRow1.siteInterval);

                            cellIds.push_back(UINT_MAX);
                            movables.push_back(cellIds);
                            cellIds.pop_back();
                            sitesPerWS = int(floor(numSitesWS/2+0.5));
                            numSitesWS -= sitesPerWS;
                            whiteSpaceWidths.push_back(sitesPerWS*
                                    smallRow1.siteInterval);
                            currWSLoc.x=smallRow1.xMax;
                            currWSLoc.y=smallRow1.yMin;
                            whiteSpaceLocs.push_back(currWSLoc);
                            currWSLoc.x += sitesPerWS*smallRow1.siteInterval ;
                            smallRow1.xMax += sitesPerWS*smallRow1.siteInterval;
                            currWSWidth -= sitesPerWS*smallRow1.siteInterval;
                            ++totalMovables1;
                            if(totalMovables1 == windowSize1)
                                break;

                            cellIds.push_back(UINT_MAX);
                            movables.push_back(cellIds);
                            cellIds.pop_back();
                            sitesPerWS = numSitesWS;
                            whiteSpaceWidths.push_back(sitesPerWS*
                                    smallRow1.siteInterval);
                            whiteSpaceLocs.push_back(currWSLoc);
                            currWSLoc.x += sitesPerWS*smallRow1.siteInterval ;
                            smallRow1.xMax += sitesPerWS*smallRow1.siteInterval;
                            currWSWidth -= sitesPerWS*smallRow1.siteInterval;
                            ++totalMovables1;
                            if(totalMovables1 == windowSize1)
                                break;
                        }
                        else if(currWSWidth > 0 && 
                                currWSWidth < 2*smallRow1.siteInterval)
                        {
                            cellIds.push_back(UINT_MAX);
                            movables.push_back(cellIds);
                            cellIds.pop_back();
                            whiteSpaceWidths.push_back(currWSWidth);
                            currWSLoc.x=smallRow1.xMax;
                            currWSLoc.y=smallRow1.yMin;
                            whiteSpaceLocs.push_back(currWSLoc);
                            currWSLoc.x += currWSWidth;
                            smallRow1.xMax += currWSWidth;
                            ++totalMovables1;
                            if(totalMovables1 == windowSize1)
                                break;
                        }
                    }
                }
                smallPbRows.push_back(smallRow1);   
                row1Size = movables.size();

                //now start the second row
                smallRow2.xMax = smallRow2.xMin = 
                    rbplace[*(subrow2.cellIdsBegin()+firstCell2)].x;

                if(rIt1 != rIt2 && windowSize2 != 0)
                {
                    for(nodeCtr2 = firstCell2; nodeCtr2 < totalCells2;)
                    {
                        vector <unsigned> cellIds;
                        //if not last cell
                        if(nodeCtr2 < totalCells2-1)
                        {
                            unsigned currentCell = *(subrow2.cellIdsBegin()+nodeCtr2);
                            unsigned nextCell = *(subrow2.cellIdsBegin()+nodeCtr2+1);
                            double totalWhiteSpace = rbplace[nextCell].x - rbplace[currentCell].x - hgWDims.getNodeWidth(currentCell);

                            //if whitespace found
                            if(totalWhiteSpace > 0)
                            {
                                cellIds.push_back(currentCell);
                                movables.push_back(cellIds);
                                cellIds.pop_back();
                                smallRow2.xMax = rbplace[currentCell].x + hgWDims.getNodeWidth(currentCell);
                                ++nodeCtr2;
                                ++totalMovables2;
                                if(totalMovables2 == windowSize2)
                                    break;

                                Point currWSLoc;
                                currWSLoc.x = smallRow2.xMax;
                                currWSLoc.y = smallRow2.yMin;
                                bool LR = true;
                                bool movablesFull = divideWhitespace(subrow2, windowSize2,
                                     LR, smallRow2, totalWhiteSpace, currWSLoc, totalMovables2, movables,
                                     whiteSpaceWidths, whiteSpaceLocs, params.wsChunkCap);

                                if(movablesFull)
                                    break;
                            }
                            else //cluster here
                            {
                                double rd = rndDBL;
                                cellIds.push_back(currentCell);
                                smallRow2.xMax = rbplace[currentCell].x + hgWDims.getNodeWidth(currentCell);
                                ++nodeCtr2;
                                if(rd > 0 && cluster == 1)
                                {
                                    cellIds.push_back(nextCell);
                                    smallRow2.xMax = rbplace[nextCell].x + hgWDims.getNodeWidth(nextCell);
                                    ++nodeCtr2;			       
                                }
                                movables.push_back(cellIds);
                                cellIds.clear();
                                ++totalMovables2;
                                if(totalMovables2 == windowSize2)
                                    break;

                                //if white space present after clustered cell
                                if(nodeCtr2 <= totalCells2-1 && rd > 0 && cluster == 1) 
                                {
                                    unsigned thirdCell = *(subrow2.cellIdsBegin()+nodeCtr2);
                                    totalWhiteSpace = rbplace[thirdCell].x - rbplace[nextCell].x - hgWDims.getNodeWidth(nextCell);

                                    Point currWSLoc;
                                    if(totalWhiteSpace >= 2*smallRow2.siteInterval)
                                    {
                                        numSitesWS = int(totalWhiteSpace/smallRow2.siteInterval);
                                        cellIds.push_back(UINT_MAX);
                                        movables.push_back(cellIds);
                                        cellIds.pop_back();
                                        sitesPerWS = int(floor(numSitesWS/2+0.5));
                                        numSitesWS -= sitesPerWS;
                                        whiteSpaceWidths.push_back(sitesPerWS*smallRow2.siteInterval);
                                        currWSLoc.x=smallRow2.xMax;
                                        currWSLoc.y=smallRow2.yMin;
                                        whiteSpaceLocs.push_back(currWSLoc);
                                        currWSLoc.x += sitesPerWS*smallRow2.siteInterval;
                                        smallRow2.xMax += sitesPerWS*smallRow2.siteInterval;
                                        totalWhiteSpace -= sitesPerWS*smallRow2.siteInterval;
                                        ++totalMovables2;
                                        if(totalMovables2 == windowSize2)
                                            break;

                                        cellIds.push_back(UINT_MAX);
                                        movables.push_back(cellIds);
                                        cellIds.pop_back();
                                        sitesPerWS = numSitesWS;
                                        whiteSpaceWidths.push_back(sitesPerWS*smallRow2.siteInterval);
                                        whiteSpaceLocs.push_back(currWSLoc);
                                        currWSLoc.x += sitesPerWS*smallRow2.siteInterval;
                                        smallRow2.xMax += sitesPerWS*smallRow2.siteInterval;
                                        totalWhiteSpace -= sitesPerWS*smallRow2.siteInterval;
                                        ++totalMovables2;
                                        if(totalMovables2 == windowSize2)
                                            break;
                                    }
                                    else if(totalWhiteSpace > 0)
                                    {
                                        cellIds.push_back(UINT_MAX);
                                        movables.push_back(cellIds);
                                        cellIds.pop_back();
                                        whiteSpaceWidths.push_back(totalWhiteSpace);
                                        currWSLoc.x=smallRow2.xMax;
                                        currWSLoc.y=smallRow2.yMin;
                                        whiteSpaceLocs.push_back(currWSLoc);
                                        currWSLoc.x += totalWhiteSpace;
                                        smallRow2.xMax += totalWhiteSpace;
                                        ++totalMovables2;
                                        if(totalMovables2 == windowSize2)
                                            break;
                                    }
                                }
                            }
                        }
                        else
                        {
                            cellIds.push_back(*(subrow2.cellIdsBegin()+nodeCtr2));
                            movables.push_back(cellIds);
                            cellIds.pop_back();
                            smallRow2.xMax=
                                rbplace[*(subrow2.cellIdsBegin()+nodeCtr2)].x+
                                hgWDims.getNodeWidth(*(subrow2.cellIdsBegin()+nodeCtr2));
                            ++nodeCtr2;
                            ++totalMovables2;
                            if(totalMovables2 == windowSize2)
                                break;		  

                            //if whitespace present after last cell
                            double currWSWidth=subrow2.getXEnd() - smallRow2.xMax;
                            Point currWSLoc;
                            if(currWSWidth >= 2*smallRow2.siteInterval)
                            {
                                numSitesWS = int(currWSWidth/smallRow2.siteInterval);

                                cellIds.push_back(UINT_MAX);
                                movables.push_back(cellIds);
                                cellIds.pop_back();
                                sitesPerWS = int(floor(numSitesWS/2+0.5));
                                numSitesWS -= sitesPerWS;
                                whiteSpaceWidths.push_back(sitesPerWS*
                                        smallRow2.siteInterval);
                                currWSLoc.x=smallRow2.xMax;
                                currWSLoc.y=smallRow2.yMin;
                                whiteSpaceLocs.push_back(currWSLoc);
                                currWSLoc.x += sitesPerWS*smallRow2.siteInterval ;
                                smallRow2.xMax += sitesPerWS*smallRow2.siteInterval;
                                currWSWidth -= sitesPerWS*smallRow2.siteInterval;
                                ++totalMovables2;
                                if(totalMovables2 == windowSize2)
                                    break;

                                cellIds.push_back(UINT_MAX);
                                movables.push_back(cellIds);
                                cellIds.pop_back();
                                sitesPerWS = numSitesWS;
                                whiteSpaceWidths.push_back(sitesPerWS*
                                        smallRow2.siteInterval);
                                whiteSpaceLocs.push_back(currWSLoc);
                                currWSLoc.x += sitesPerWS*smallRow2.siteInterval ;
                                smallRow2.xMax += sitesPerWS*smallRow2.siteInterval;
                                currWSWidth -= sitesPerWS*smallRow2.siteInterval;
                                ++totalMovables2;
                                if(totalMovables2 == windowSize2)
                                    break;
                            }
                            else if(currWSWidth > 0 && 
                                    currWSWidth < 2*smallRow2.siteInterval)
                            {
                                cellIds.push_back(UINT_MAX);
                                movables.push_back(cellIds);
                                cellIds.pop_back();
                                whiteSpaceWidths.push_back(currWSWidth);
                                currWSLoc.x=smallRow2.xMax;
                                currWSLoc.y=smallRow2.yMin;
                                whiteSpaceLocs.push_back(currWSLoc);
                                currWSLoc.x += currWSWidth;
                                smallRow2.xMax += currWSWidth;
                                ++totalMovables2;
                                if(totalMovables2 == windowSize2)
                                    break;
                            }
                        }
                    }
                    smallPbRows.push_back(smallRow2);
                }

                if(movables.size() < 2)
                {
                    firstCell1=nodeCtr1;
                    continue;
                }

                currCellLoc1=nodeCtr1;
                currCellLoc2=nodeCtr2;

                //assume all the cells are placed
                //assume all cells are not fixed

                bit_vector fixed(movables.size(),0);
                for(unsigned i=0; i<movables.size(); ++i)
                {
                    bool oneFixedCell=false;
                    for(unsigned j=0; j<movables[i].size(); ++j)
                    {
                        if(movables[i][j] != UINT_MAX && rbplace.isFixed(movables[i][j]))
                        {
                            oneFixedCell = true;
                            break;
                        }
                    }
                    if(oneFixedCell)
                        fixed[i] = true;
                }

                // the new placement of the movables will go here
                Placement soln(movables.size());

#ifdef USEFLUTE
                if(params.useSteiner && smallPbRows.size() == 1)
                {
                  if(params.useDP)
                  {
                    placeWithFluteDP(rbplace, movables, smallPbRows[0],
                                     whiteSpaceWidths, forFlute, soln, tempSolnForFlute, ptSets);
                  }
                  else
                  {
                    placeWithFlute(rbplace, movables, smallPbRows[0],
                                   whiteSpaceWidths, forFlute, soln, tempSolnForFlute, ptSets);
                  }
                }
                else
#endif
                { // royj: generates a solution for this set of movables using a SmallPlacer
                  SmallPlacementProblem smallProblem(hgWDims, rbplace,
                          placed, fixed, movables, smallPbRows, whiteSpaceWidths, spBitBoards);

                  SmallPlacementSolution& solution=smallProblem.getSolnBuffer()[0];
                  solution.populated=true;

                  unsigned currWSLocation=0;
                  for(unsigned cell=0; cell!=movables.size(); cell++)
                  {
                      if(movables[cell][0] != UINT_MAX)
                      {
                          solution.placement[cell]   =rbplace[movables[cell][0]];
                          solution.orientations[cell]=rbplace.getOrient(movables[cell][0]);
                      }
                      else
                      {
                          Point WSLoc = whiteSpaceLocs[currWSLocation];
                          ++currWSLocation;
                          solution.placement[cell] = WSLoc;
                      }
                  }

                  if(smallProblem.getNumNets() > 63)
                  {
                      ++largeNetProblems;
                      if(largeNetProblems < 4) abkwarn(0,"skipping problem with > 63 nets\n");
                      if(largeNetProblems == 3) abkwarn(0,"all further > 63 nets messages will be supressed\n");
                      firstCell1 = nodeCtr1;
                      continue;
                  }

                  //invoke the SmallPlacer to solve the problem
                  //smallProblem.plot("smPbIn");
                  //smallProblem.save("smPbIn");
                  SmallPlacer placer(smallProblem,params.smplParams,maxMem);
                  //smallProblem.plot("smPbOut");
                  //smallProblem.save("smPbOut");

                  //get best solution
                  //soln = smallProblem.getBestSoln().placement;
                  abkassert(smallProblem.getBestSoln().placement.getSize() == soln.size(),"undersized placement");

                  for(unsigned i = 0; i < soln.size(); ++i)
                  {
                    soln[i] = smallProblem.getBestSoln().placement[i];
                  }
                } // royj: end of solution generation by SmallPlacer


                //rbplace.saveAsPlot("plotWNames",-DBL_MAX,DBL_MAX,-DBL_MAX,DBL_MAX,"RBPlIn");
                //double initHPWL = rbplace.evalHPWL(true);
                rbplace.updateIronedCellsLR(movables,soln,
                        const_cast<RBPSubRow&>(subrow1),
                        const_cast<RBPSubRow&>(subrow2) );

#ifdef USEFLUTE
                // update tempSolnForFlute so it matches the placement in rbplace
                for(unsigned i = 0; i < movables.size(); ++i)
                {
                  for(unsigned j = 0; j < movables[i].size(); ++j)
                  {
                    if(movables[i][j] != UINT_MAX)
                    {
                      const unsigned cellId = movables[i][j];
                      const Point &oldPos = tempSolnForFlute[cellId];
                      const Point &newPos = rbplace.getPlacement()[cellId];

                      if(oldPos != newPos)
                      {
                        const HGFNode &cell = hgWDims.getNodeByIdx(cellId);
                        unsigned offset = 0;
                        for(itHGFEdgeLocal e = cell.edgesBegin(); e != cell.edgesEnd(); ++e, ++offset)
                        {
                          unsigned netId = (*e)->getIndex();
                          BBox pinOffset = hgWDims.edgesOnNodePinOffsetBBox(offset, cellId, rbplace.getOrient(cellId));

                          removeFromPointSet(ptSets[netId], oldPos + pinOffset.getBottomLeft());
                          addToPointSet(ptSets[netId], newPos + pinOffset.getBottomLeft());

                          if(!equalDouble(pinOffset.xMin,pinOffset.xMax) ||
                             !equalDouble(pinOffset.yMin,pinOffset.yMax))
                          {
                            removeFromPointSet(ptSets[netId], oldPos + pinOffset.getTopRight());
                            addToPointSet(ptSets[netId], newPos + pinOffset.getTopRight());
                          }
                        }

                        tempSolnForFlute[cellId] = newPos;
                      }
                    }
                  }
                }
#endif

                // now that the cell locations have been changed in the subrows,
                // the subrow cellid lists have most likely changed
                // (potentially different order and size)
                // take care not to use outdated information

                //double finalHPWL = rbplace.evalHPWL(true);
                /*
                   if(finalHPWL > initHPWL)
                   {
                   rbplace.saveAsPlot("plotWNames",-DBL_MAX,DBL_MAX,-DBL_MAX,DBL_MAX,"RBPlOut");
                   cout<<"increase in HPWL is "<<finalHPWL-initHPWL<<endl;
                   abkfatal(0,"bye\n");
                   }
                */
                nodeCtr1=firstCell1;
                for(unsigned i=0;i<movables.size();++i)
                {
                    if(movables[i][0] != UINT_MAX)
                        if(soln[i].y == subrow1.getYCoord())
                            ++firstCell1;
                }

                //to avoid a deadlock condition
                unsigned overlap=0;
                if(row1Size == 0)
                    overlap = 0;
                else if(params.overlap > row1Size)
                    overlap = row1Size-1;
                else
                    overlap = params.overlap;

                if(nodeCtr1 < (firstCell1 - overlap))
                    firstCell1 = (firstCell1 - overlap);

                totalCells1 = subrow1.getNumCells();
                totalCells2 = subrow2.getNumCells();
            }

            if(totalCells1 < 1 || totalCells2 < 1) continue;      
        }

        if(verticalDir == 0)
        {
            ++rIt1;
            if(rIt1 == rbplace.coreRowsEnd())
                break;
        }
        else
        {
            if(rIt1 == rbplace.coreRowsBegin())
                break;
            else
                --rIt1;
        }
    }

    if(twoDim) {
        params.smplParams.algo = algo_save;
    }
}


//to iron Rows in two dimension right to left
void ironRows2DRL(RBPlacement& rbplace, const RowIroningParameters& params, SmallPlacementBitBoardContainer &spBitBoards,
		  MaxMem *maxMem, bool verticalDir, bool twoDim, bool cluster)
{
    const HGraphWDimensions& hgWDims = rbplace.getHGraph();
#ifdef USEFLUTE
    BitBoard forFlute(hgWDims.getNumEdges());
    Placement tempSolnForFlute(rbplace.getPlacement());

    vector< vector<PtSetPoint> > ptSets(hgWDims.getNumEdges());

    // populate the point sets
    for(unsigned i = 0; i < ptSets.size(); ++i)
    {
      const HGFEdge &net = hgWDims.getEdgeByIdx(i);
      ptSets[i].reserve(2*net.getDegree());

      unsigned offset = 0;
      for(itHGFNodeLocal c = net.nodesBegin(); c != net.nodesEnd(); ++c, ++offset)
      {
        unsigned cellId = (*c)->getIndex();
        BBox pinOffset = hgWDims.nodesOnEdgePinOffsetBBox(offset, i, rbplace.getOrient(cellId));

        addToPointSet(ptSets[i], tempSolnForFlute[cellId] + pinOffset.getBottomLeft());

        if(!equalDouble(pinOffset.xMin,pinOffset.xMax) ||
           !equalDouble(pinOffset.yMin,pinOffset.yMax))
        {
          addToPointSet(ptSets[i], tempSolnForFlute[cellId] + pinOffset.getTopRight());
        }
      }
    }
#endif
    bit_vector placed = bit_vector(rbplace.getNumCells(),1);

//    char fileName1[15];
//    char fileName2[15];
//    char fileName3[15];
//    char fileName4[15];
//    unsigned instance = 0;
    unsigned largeNetProblems = 0;

    RandomDouble rndDBL(-1,1);

    SmallPlParams::SmPlAlgo algo_save;
    if(twoDim) {
        algo_save = params.smplParams.algo;
        params.smplParams.algo = SmallPlParams::DynamicProgram;
    }

    //verticalDir=0 means bottom up pass else
    //top down pass
    itRBPCoreRow rIt1,rIt2;
    itRBPSubRow sIt1,sIt2;

    if(verticalDir == 0)
        rIt1=rbplace.coreRowsBegin();
    else
        rIt1=rbplace.coreRowsEnd()-1;

    while(true)                 //exit only by a break
    {
        unsigned windowSize1=params.windowSize;
        unsigned windowSize2=0;
        rIt2 = rIt1;
        if(verticalDir == 0)
        {
            if(rIt1 != rbplace.coreRowsEnd()-1  && twoDim == 1) 
                //if not last core row
            {
                ++rIt2;
                windowSize1 = unsigned(ceil(params.windowSize/2.0));
                windowSize2 = params.windowSize - windowSize1;
            }
        }
        else
        {
            if(rIt1 != rbplace.coreRowsBegin() && twoDim == 1)   
                //if not first Row
            {
                --rIt2;
                windowSize1 = unsigned(ceil(params.windowSize/2.0));
                windowSize2 = params.windowSize - windowSize1;
            }
        }

        sIt1=rIt1->subRowsEnd()-1;
        sIt2=rIt2->subRowsEnd()-1;

        while(true)
        {
            const RBPSubRow& subrow1=*sIt1;
            const RBPSubRow& subrow2=*sIt2;

            SmallPlacementRow smallRow1;
            SmallPlacementRow smallRow2;

            smallRow1.siteInterval=subrow1.getSpacing();
            smallRow1.xMax =smallRow1.xMin=subrow1.getXEnd();
            smallRow1.yMax =smallRow1.yMin=subrow1.getYCoord();
            smallRow1.yMax+=subrow1.getHeight();
            unsigned firstCell1=0;
            unsigned totalCells1=subrow1.getNumCells();
            unsigned currCellLoc1=0;
            int nodeCtr1=0;
            int lastCell1 = totalCells1-1;

            smallRow2.siteInterval=subrow2.getSpacing();
            smallRow2.xMax =smallRow2.xMin=subrow2.getXEnd();
            smallRow2.yMax =smallRow2.yMin=subrow2.getYCoord();
            smallRow2.yMax+=subrow2.getHeight();
            //       unsigned firstCell2=0;
            unsigned totalCells2=subrow2.getNumCells();
            unsigned currCellLoc2=0;
            int nodeCtr2=0;
            int lastCell2 = totalCells2-1;

            int numSitesWS=0;
            int sitesPerWS=0;

            if(totalCells1 < 1 || totalCells2 < 1) 
            {
                if(sIt1 == rIt1->subRowsBegin() || sIt2 == rIt2->subRowsBegin())
                    break;
                --sIt1;
                --sIt2;
                continue;
            }

            //now start forming instances
            while (lastCell1 >= 0)
            { 
//                ++instance;
//                sprintf(fileName1,"SmPbIn-%d",instance);
//                sprintf(fileName2,"SmPbOut-%d",instance);
//                sprintf(fileName3,"RbPbIn-%d",instance);
//                sprintf(fileName4,"RbPbOut-%d",instance);

                vector<SmallPlacementRow> smallPbRows;
                //form the 1'st row first
                //get the cells + whitespace cells here
                vector< vector<unsigned> > movables;
                vector<double> whiteSpaceWidths;
                vector<Point> whiteSpaceLocs;
                unsigned totalMovables1=0;
                unsigned totalMovables2=0;
                unsigned row1Size=0;

                bool brokeFromLoop = 0;
                if(lastCell1 >= 0 )
                {
                    for(nodeCtr2=0;nodeCtr2<int(totalCells2-1);++nodeCtr2)
                    {
                        if((rbplace[*(subrow2.cellIdsBegin()+nodeCtr2)].x <=
                                    rbplace[*(subrow1.cellIdsBegin()+lastCell1)].x) &&
                                (rbplace[*(subrow2.cellIdsBegin()+nodeCtr2+1)].x >
                                 rbplace[*(subrow1.cellIdsBegin()+lastCell1)].x))
                        {
                            brokeFromLoop = 1;
                            break;
                        }
                    }
                    if(!brokeFromLoop)
                        nodeCtr2 = totalCells2-1;
                    lastCell2 = nodeCtr2;
                }
                else
                    lastCell2 = -1;

                nodeCtr1=lastCell1;
                nodeCtr2=lastCell2;

                smallRow1.xMax = smallRow1.xMin = 
                    rbplace[*(subrow1.cellIdsBegin()+lastCell1)].x+
                    hgWDims.getNodeWidth(*(subrow1.cellIdsBegin()+lastCell1));

                while(nodeCtr1 >=0)
                {
                    vector<unsigned> cellIds;

                    //if not last cell
                    if(nodeCtr1 > 0)
                    {
                        unsigned currentCell = *(subrow1.cellIdsBegin()+nodeCtr1);
                        unsigned prevCell = *(subrow1.cellIdsBegin()+nodeCtr1-1);
                        double totalWhiteSpace = rbplace[currentCell].x - rbplace[prevCell].x - hgWDims.getNodeWidth(prevCell);

                        //if whitespace found
                        if(totalWhiteSpace > 0)
                        {
                            cellIds.push_back(currentCell);
                            movables.push_back(cellIds);
                            cellIds.pop_back();
                            smallRow1.xMin = rbplace[currentCell].x;
                            --nodeCtr1;
                            ++totalMovables1;
                            if(totalMovables1 == windowSize1)
                                break;

                            Point currWSLoc;
                            currWSLoc.x = smallRow1.xMin;
                            currWSLoc.y = smallRow1.yMin;
                            bool LR = false;
                            bool movablesFull = divideWhitespace(subrow1, windowSize1,
                                 LR, smallRow1, totalWhiteSpace, currWSLoc, totalMovables1, movables,
                                 whiteSpaceWidths, whiteSpaceLocs, params.wsChunkCap);

                            if(movablesFull)
                                break;
                        }
                        else  //cluster the cells here
                        {
                            double rd = rndDBL;
                            cellIds.push_back(currentCell);
                            smallRow1.xMin = rbplace[currentCell].x;
                            --nodeCtr1;
                            //cluster the next cell based on random value
                            if(rd > 0 && cluster == 1)  //cluster
                            {
                                //cellIds.push_back(*(subrow1.cellIdsBegin()+nodeCtr1));
                                cellIds.insert(cellIds.begin(), prevCell);
                                smallRow1.xMin = rbplace[prevCell].x;
                                --nodeCtr1; 
                            }

                            movables.push_back(cellIds);
                            cellIds.clear();
                            ++totalMovables1;
                            if(totalMovables1 == windowSize1)
                                break;

                            //if white space present after clustered cell
                            if(nodeCtr1 >= 0 && rd > 0 && cluster == 1)
                                //if not first cell
                            {
                                unsigned thirdCell = *(subrow1.cellIdsBegin()+nodeCtr1);
                                totalWhiteSpace = rbplace[prevCell].x - rbplace[thirdCell].x - hgWDims.getNodeWidth(thirdCell);

                                Point currWSLoc;
                                if(totalWhiteSpace >= 2*smallRow1.siteInterval)
                                {
                                    numSitesWS = int(totalWhiteSpace/smallRow1.siteInterval);

                                    cellIds.push_back(UINT_MAX);
                                    movables.push_back(cellIds);
                                    cellIds.pop_back();
                                    sitesPerWS = int(floor(numSitesWS/2+0.5));
                                    numSitesWS -= sitesPerWS;
                                    whiteSpaceWidths.push_back(sitesPerWS*smallRow1.siteInterval);
                                    currWSLoc.x=smallRow1.xMin;
                                    currWSLoc.y=smallRow1.yMin;
                                    currWSLoc.x -= sitesPerWS*smallRow1.siteInterval;
                                    whiteSpaceLocs.push_back(currWSLoc);
                                    smallRow1.xMin -= sitesPerWS*smallRow1.siteInterval;
                                    totalWhiteSpace -= sitesPerWS*smallRow1.siteInterval;
                                    ++totalMovables1;
                                    if(totalMovables1 == windowSize1)
                                        break;

                                    cellIds.push_back(UINT_MAX);
                                    movables.push_back(cellIds);
                                    cellIds.pop_back();
                                    sitesPerWS = numSitesWS;
                                    whiteSpaceWidths.push_back(sitesPerWS*smallRow1.siteInterval);
                                    currWSLoc.x -= sitesPerWS*smallRow1.siteInterval;
                                    whiteSpaceLocs.push_back(currWSLoc);
                                    smallRow1.xMin -= sitesPerWS*smallRow1.siteInterval;
                                    totalWhiteSpace -= sitesPerWS*smallRow1.siteInterval;
                                    ++totalMovables1;
                                    if(totalMovables1 == windowSize1)
                                        break;
                                }
                                else if(totalWhiteSpace > 0)
                                {
                                    cellIds.push_back(UINT_MAX);
                                    movables.push_back(cellIds);
                                    cellIds.pop_back();
                                    whiteSpaceWidths.push_back(totalWhiteSpace);
                                    currWSLoc.x=smallRow1.xMin;
                                    currWSLoc.y=smallRow1.yMin;
                                    currWSLoc.x -= totalWhiteSpace;
                                    whiteSpaceLocs.push_back(currWSLoc);
                                    smallRow1.xMin -= totalWhiteSpace;
                                    ++totalMovables1;
                                    if(totalMovables1 == windowSize1)
                                        break;
                                }
                            }
                        }
                    }
                    else
                    {
                        cellIds.push_back(*(subrow1.cellIdsBegin()+nodeCtr1));
                        movables.push_back(cellIds);
                        cellIds.pop_back();
                        smallRow1.xMin = rbplace[*(subrow1.cellIdsBegin()+nodeCtr1)].x;
                        --nodeCtr1;
                        ++totalMovables1;
                        if(totalMovables1 == windowSize1)
                            break;		  

                        //if whitespace present after last cell
                        double currWSWidth=subrow1.getXStart() - smallRow1.xMin;
                        Point currWSLoc;
                        if(currWSWidth >= 2*smallRow1.siteInterval)
                        {
                            numSitesWS = int(currWSWidth/smallRow1.siteInterval);

                            cellIds.push_back(UINT_MAX);
                            movables.push_back(cellIds);
                            cellIds.pop_back();
                            sitesPerWS = int(floor(numSitesWS/2+0.5));
                            numSitesWS -= sitesPerWS;
                            whiteSpaceWidths.push_back(sitesPerWS*
                                    smallRow1.siteInterval);
                            currWSLoc.x=smallRow1.xMin;
                            currWSLoc.y=smallRow1.yMin;
                            currWSLoc.x -= sitesPerWS*smallRow1.siteInterval ;
                            whiteSpaceLocs.push_back(currWSLoc);
                            smallRow1.xMin -= sitesPerWS*smallRow1.siteInterval;
                            currWSWidth -= sitesPerWS*smallRow1.siteInterval;
                            ++totalMovables1;
                            if(totalMovables1 == windowSize1)
                                break;

                            cellIds.push_back(UINT_MAX);
                            movables.push_back(cellIds);
                            cellIds.pop_back();
                            sitesPerWS = numSitesWS;
                            whiteSpaceWidths.push_back(sitesPerWS*
                                    smallRow1.siteInterval);
                            currWSLoc.x -= sitesPerWS*smallRow1.siteInterval ;
                            whiteSpaceLocs.push_back(currWSLoc);
                            smallRow1.xMin -= sitesPerWS*smallRow1.siteInterval;
                            currWSWidth -= sitesPerWS*smallRow1.siteInterval;
                            ++totalMovables1;
                            if(totalMovables1 == windowSize1)
                                break;
                        }
                        else if(currWSWidth > 0 && 
                                currWSWidth < 2*smallRow1.siteInterval)
                        {
                            cellIds.push_back(UINT_MAX);
                            movables.push_back(cellIds);
                            cellIds.pop_back();
                            whiteSpaceWidths.push_back(currWSWidth);
                            currWSLoc.x=smallRow1.xMin;
                            currWSLoc.y=smallRow1.yMin;
                            currWSLoc.x -= currWSWidth;
                            whiteSpaceLocs.push_back(currWSLoc);
                            smallRow1.xMin -= currWSWidth;
                            ++totalMovables1;
                            if(totalMovables1 == windowSize1)
                                break;
                        }
                    }
                }
                smallPbRows.push_back(smallRow1);   
                row1Size = movables.size();
                firstCell1 = nodeCtr1+1;

                //now start the second row
                nodeCtr2=lastCell2;
                smallRow2.xMax = smallRow2.xMin = 
                    rbplace[*(subrow2.cellIdsBegin()+lastCell2)].x+
                    hgWDims.getNodeWidth(*(subrow2.cellIdsBegin()+lastCell2));

                if(rIt1 != rIt2 && windowSize2 != 0)
                {
                    while(nodeCtr2 >= 0)
                    {
                        vector<unsigned> cellIds;
                        //if not last cell
                        if(nodeCtr2 > 0)
                        {
                            unsigned currentCell = *(subrow2.cellIdsBegin()+nodeCtr2);
                            unsigned prevCell = *(subrow2.cellIdsBegin()+nodeCtr2-1);
                            double totalWhiteSpace = rbplace[currentCell].x - rbplace[prevCell].x - hgWDims.getNodeWidth(prevCell);

                            //if whitespace found
                            if(totalWhiteSpace > 0)
                            {
                                cellIds.push_back(currentCell);
                                movables.push_back(cellIds);
                                cellIds.pop_back();
                                smallRow2.xMin = rbplace[currentCell].x;
                                --nodeCtr2;
                                ++totalMovables2;
                                if(totalMovables2 == windowSize2)
                                    break;

                                Point currWSLoc;
                                currWSLoc.x = smallRow2.xMin;
                                currWSLoc.y = smallRow2.yMin;
                                bool LR = false;
                                bool movablesFull = divideWhitespace(subrow2, windowSize2,
                                     LR, smallRow2, totalWhiteSpace, currWSLoc, totalMovables2, movables,
                                     whiteSpaceWidths, whiteSpaceLocs, params.wsChunkCap);

                                if(movablesFull)
                                    break;  
                            }
                            else    //cluster here
                            {
                                double rd = rndDBL;
                                cellIds.push_back(currentCell);
                                smallRow2.xMin = rbplace[currentCell].x;
                                --nodeCtr2;
                                //cluster the next cell based on random value
                                if(rd > 0 && cluster == 1)  //cluster
                                {
                                    //cellIds.push_back(prevCell);
                                    cellIds.insert(cellIds.begin(),prevCell);

                                    smallRow2.xMin = rbplace[prevCell].x;
                                    --nodeCtr2; 
                                }

                                movables.push_back(cellIds);
                                cellIds.clear();
                                ++totalMovables2;
                                if(totalMovables2 == windowSize2)
                                    break;

                                //if white space present after clustered cell
                                if(nodeCtr2 >= 0 && rd > 0 && cluster == 1)
                                    //if not beyond first cell
                                {
                                    unsigned thirdCell = *(subrow2.cellIdsBegin()+nodeCtr2);
                                    totalWhiteSpace = rbplace[prevCell].x - rbplace[thirdCell].x - hgWDims.getNodeWidth(thirdCell);

                                    Point currWSLoc;
                                    if(totalWhiteSpace >= 2*smallRow2.siteInterval)
                                    {
                                        numSitesWS = int(totalWhiteSpace/smallRow2.siteInterval);
                                        cellIds.push_back(UINT_MAX);
                                        movables.push_back(cellIds);
                                        cellIds.pop_back();
                                        sitesPerWS = int(floor(numSitesWS/2+0.5));
                                        numSitesWS -= sitesPerWS;
                                        whiteSpaceWidths.push_back(sitesPerWS*smallRow2.siteInterval);
                                        currWSLoc.x=smallRow2.xMin;
                                        currWSLoc.y=smallRow2.yMin;
                                        currWSLoc.x -= sitesPerWS*smallRow2.siteInterval ;
                                        whiteSpaceLocs.push_back(currWSLoc);
                                        smallRow2.xMin -= sitesPerWS*smallRow2.siteInterval;
                                        totalWhiteSpace -= sitesPerWS*smallRow2.siteInterval;
                                        ++totalMovables2;
                                        if(totalMovables2 == windowSize2)
                                            break;

                                        cellIds.push_back(UINT_MAX);
                                        movables.push_back(cellIds);
                                        cellIds.pop_back();
                                        sitesPerWS = numSitesWS;
                                        whiteSpaceWidths.push_back(sitesPerWS*smallRow2.siteInterval);
                                        currWSLoc.x -= sitesPerWS*smallRow2.siteInterval;
                                        whiteSpaceLocs.push_back(currWSLoc);
                                        smallRow2.xMin -= sitesPerWS*smallRow2.siteInterval;
                                        totalWhiteSpace -= sitesPerWS*smallRow2.siteInterval;
                                        ++totalMovables2;
                                        if(totalMovables2 == windowSize2)
                                            break;
                                    }
                                    else if(totalWhiteSpace > 0)
                                    {
                                        cellIds.push_back(UINT_MAX);
                                        movables.push_back(cellIds);
                                        cellIds.pop_back();
                                        whiteSpaceWidths.push_back(totalWhiteSpace);
                                        currWSLoc.x=smallRow2.xMin;
                                        currWSLoc.y=smallRow2.yMin;
                                        currWSLoc.x -= totalWhiteSpace;
                                        whiteSpaceLocs.push_back(currWSLoc);
                                        smallRow2.xMin -= totalWhiteSpace;
                                        ++totalMovables2;
                                        if(totalMovables2 == windowSize2)
                                            break;
                                    }
                                }
                            }
                        }
                        else
                        {
                            cellIds.push_back(*(subrow2.cellIdsBegin()+nodeCtr2));
                            movables.push_back(cellIds);
                            cellIds.pop_back();
                            smallRow2.xMin=
                                rbplace[*(subrow2.cellIdsBegin()+nodeCtr2)].x;
                            --nodeCtr2;
                            ++totalMovables2;
                            if(totalMovables2 == windowSize2)
                                break;		  

                            //if whitespace present after last cell
                            double currWSWidth=subrow2.getXStart() - smallRow2.xMin;
                            Point currWSLoc;
                            if(currWSWidth >= 2*smallRow2.siteInterval)
                            {
                                numSitesWS = int(currWSWidth/smallRow2.siteInterval);
                                cellIds.push_back(UINT_MAX);
                                movables.push_back(cellIds);
                                cellIds.pop_back();
                                sitesPerWS = int(floor(numSitesWS/2+0.5));
                                numSitesWS -= sitesPerWS;
                                whiteSpaceWidths.push_back(sitesPerWS*
                                        smallRow2.siteInterval);
                                currWSLoc.x=smallRow2.xMin;
                                currWSLoc.y=smallRow2.yMin;
                                currWSLoc.x -= sitesPerWS*smallRow2.siteInterval ;
                                whiteSpaceLocs.push_back(currWSLoc);
                                smallRow2.xMin -= sitesPerWS*smallRow2.siteInterval;
                                currWSWidth -= sitesPerWS*smallRow2.siteInterval;
                                ++totalMovables2;
                                if(totalMovables2 == windowSize2)
                                    break;

                                cellIds.push_back(UINT_MAX);
                                movables.push_back(cellIds);
                                cellIds.pop_back();
                                sitesPerWS = numSitesWS;
                                whiteSpaceWidths.push_back(sitesPerWS*
                                        smallRow2.siteInterval);
                                currWSLoc.x -= sitesPerWS*smallRow2.siteInterval ;
                                whiteSpaceLocs.push_back(currWSLoc);
                                smallRow2.xMin -= sitesPerWS*smallRow2.siteInterval;
                                currWSWidth -= sitesPerWS*smallRow2.siteInterval;
                                ++totalMovables2;
                                if(totalMovables2 == windowSize2)
                                    break;
                            }
                            else if(currWSWidth > 0 && 
                                    currWSWidth < 2*smallRow2.siteInterval)
                            {
                                cellIds.push_back(UINT_MAX);
                                movables.push_back(cellIds);
                                cellIds.pop_back();
                                whiteSpaceWidths.push_back(currWSWidth);
                                currWSLoc.x=smallRow2.xMin;
                                currWSLoc.y=smallRow2.yMin;
                                currWSLoc.x -= currWSWidth;
                                whiteSpaceLocs.push_back(currWSLoc);
                                smallRow2.xMin -= currWSWidth;
                                ++totalMovables2;
                                if(totalMovables2 == windowSize2)
                                    break;
                            }
                        }
                    }
                    smallPbRows.push_back(smallRow2);
                }

                if(movables.size() < 2)
                {
                    if(nodeCtr1 == -1)
                        break;
                    lastCell1=nodeCtr1;
                    continue;
                }

                currCellLoc1=nodeCtr1;
                currCellLoc2=nodeCtr2;

                //assume all the cells are placed
                //assume all cells are not fixed
                vector< vector<unsigned> > movablesLR;

                vector<double> whiteSpaceWidthsLR;
                vector<Point> whiteSpaceLocsLR;
                for(unsigned i=1;i<=movables.size();++i)
                {
                    movablesLR.push_back(movables[movables.size()-i]);
                    /*
                       if(movables[i-1][0] != UINT_MAX)
                       cout<<hgWDims.getNodeByIdx(movables[i-1][0])<<endl;
                       else
                       cout<<UINT_MAX<<endl;
                     */
                }
                for(unsigned i=1;i<=whiteSpaceWidths.size();++i)
                {
                    whiteSpaceWidthsLR.push_back(whiteSpaceWidths
                            [whiteSpaceWidths.size()-i]);
                }
                for(unsigned i=1;i<=whiteSpaceLocs.size();++i)
                {
                    whiteSpaceLocsLR.push_back(whiteSpaceLocs
                            [whiteSpaceLocs.size()-i]);
                }

                bit_vector fixed(movablesLR.size(),0);
                for(unsigned i=0; i<movablesLR.size(); ++i)
                {
                    bool oneFixedCell=false;
                    for(unsigned j=0; j<movablesLR[i].size(); ++j)
                    {

                        if(movablesLR[i][j]!=UINT_MAX && rbplace.isFixed(movablesLR[i][j]))
                        {
                            oneFixedCell = true;
                            break;
                        }
                    }
                    if(oneFixedCell)
                        fixed[i] = true;
                }

                //get the no: of cells in initial instance
                unsigned inst1PrevSize=0;
                for(unsigned cell=0;cell<movablesLR.size();++cell)
                {
                    if(movablesLR[cell][0] != UINT_MAX)
                    {
                       const Point &posOfCell = rbplace[movablesLR[cell].back()];
                       if(posOfCell.y == subrow1.getYCoord())
                       {
                           ++inst1PrevSize;
                       }
                    }
                }

                Placement soln(movablesLR.size());

#ifdef USEFLUTE
                if(params.useSteiner && smallPbRows.size() == 1)
                {
                  if(params.useDP)
                  {
                    placeWithFluteDP(rbplace, movablesLR, smallPbRows[0],
                                     whiteSpaceWidthsLR, forFlute, soln, tempSolnForFlute, ptSets);
                  }
                  else
                  {
                    placeWithFlute(rbplace, movablesLR, smallPbRows[0],
                                   whiteSpaceWidthsLR, forFlute, soln, tempSolnForFlute, ptSets);
                  }
                }
                else
#endif
                { // royj: generates a solution for this set of movables using a SmallPlacer
                  //reverse the smallRows order to confirm with what SmallPlacer expects
                  std::reverse(smallPbRows.begin(),smallPbRows.end());

                  SmallPlacementProblem smallProblem(hgWDims, rbplace,
                          placed, fixed, movablesLR, smallPbRows, whiteSpaceWidthsLR, spBitBoards);

                  SmallPlacementSolution& solution=smallProblem.getSolnBuffer()[0];
                  solution.populated=true;

                  unsigned currWSLocation=0;
                  for(unsigned cell=0; cell!=movablesLR.size(); cell++)
                  {
                      if(movablesLR[cell][0] != UINT_MAX)
                      {
                          solution.placement[cell]   =rbplace[movablesLR[cell].back()];
                          solution.orientations[cell]=rbplace.getOrient(movablesLR[cell].back());
                      }
                      else
                      {
                          Point WSLoc = whiteSpaceLocsLR[currWSLocation];
                          ++currWSLocation;
                          solution.placement[cell] = WSLoc;
                      }
                  }

                  if(smallProblem.getNumNets() > 63)
                  {
                      ++largeNetProblems;
                      if(largeNetProblems < 4) abkwarn(0,"skipping problem with > 63 nets\n");
                      if(largeNetProblems == 3) abkwarn(0,"all further > 63 nets messages will be supressed\n");
                      lastCell1 = nodeCtr1;
                      continue;
                  }

                  //invoke the SmallPlacer to solve the problem
                  //smallProblem.plot("SmPbIn");
                  //smallProblem.save("SmPbIn");
                  SmallPlacer placer(smallProblem,params.smplParams,maxMem);
                  //smallProblem.plot("SmPbOut");
                  //smallProblem.save("SmPbOut");

                  //get best solution
                  //const Placement& soln = smallProblem.getBestSoln().placement;
                  abkassert(smallProblem.getBestSoln().placement.getSize() == soln.size(),"undersized placement");

                  for(unsigned i = 0; i < soln.size(); ++i)
                  {
                    soln[i] = smallProblem.getBestSoln().placement[i];
                  }
                } // royj: end of solution generation by SmallPlacer
                

                //rbplace.saveAsPlot("plotWNames",-DBL_MAX,DBL_MAX,-DBL_MAX,DBL_MAX,fileName3);
                //double initHPWL = rbplace.evalHPWL(true);

                rbplace.updateIronedCellsRL(movablesLR,soln,
                        const_cast<RBPSubRow&>(subrow1),
                        const_cast<RBPSubRow&>(subrow2) );

#ifdef USEFLUTE
                // update tempSolnForFlute so it matches the placement in rbplace
                for(unsigned i = 0; i < movables.size(); ++i)
                {
                  for(unsigned j = 0; j < movables[i].size(); ++j)
                  {
                    if(movables[i][j] != UINT_MAX)
                    {
                      const unsigned cellId = movables[i][j];
                      const Point &oldPos = tempSolnForFlute[cellId];
                      const Point &newPos = rbplace.getPlacement()[cellId];

                      if(oldPos != newPos)
                      {
                        const HGFNode &cell = hgWDims.getNodeByIdx(cellId);
                        unsigned offset = 0;
                        for(itHGFEdgeLocal e = cell.edgesBegin(); e != cell.edgesEnd(); ++e, ++offset)
                        {
                          unsigned netId = (*e)->getIndex();
                          BBox pinOffset = hgWDims.edgesOnNodePinOffsetBBox(offset, cellId, rbplace.getOrient(cellId));

                          removeFromPointSet(ptSets[netId], oldPos + pinOffset.getBottomLeft());
                          addToPointSet(ptSets[netId], newPos + pinOffset.getBottomLeft());

                          if(!equalDouble(pinOffset.xMin,pinOffset.xMax) ||
                             !equalDouble(pinOffset.yMin,pinOffset.yMax))
                          {
                            removeFromPointSet(ptSets[netId], oldPos + pinOffset.getTopRight());
                            addToPointSet(ptSets[netId], newPos + pinOffset.getTopRight());
                          }
                        }

                        tempSolnForFlute[cellId] = newPos;
                      }
                    }
                  }
                }
#endif

                //double finalHPWL = rbplace.evalHPWL(true);
                /*
                   if(finalHPWL > initHPWL)
                   {
                   cout<<"increase in HPWL is "<<finalHPWL-initHPWL<<endl;
                   rbplace.saveAsPlot("plotWNames",-DBL_MAX,DBL_MAX,-DBL_MAX,DBL_MAX,"RBPlOut");
                   abkfatal(0,"bye\n");
                   }
                 */
                nodeCtr1=lastCell1;
                //get the no: of cells in the final instance
                lastCell1 = lastCell1 - inst1PrevSize;
                unsigned inst1Size=0;
                for(unsigned i=0;i<movablesLR.size();++i)
                {
                    if(movablesLR[i][0] != UINT_MAX)
                        if(soln[i].y == subrow1.getYCoord())
                        {
                            ++inst1Size;
                            --lastCell1;
                        }
                }
                lastCell1 = lastCell1 + inst1Size;
                nodeCtr1=lastCell1;

                //to avoid a deadlock condition

                unsigned overlap=0;
                if(inst1Size == 0)
                    overlap = 0;
                else if(params.overlap >= inst1Size)
                    overlap = inst1Size-1;
                else
                    overlap = params.overlap;
                if(nodeCtr1 > int(lastCell1 + overlap))
                    lastCell1 = (lastCell1 + overlap);

                totalCells1=subrow1.getNumCells();
                totalCells2=subrow2.getNumCells();
            }

            if(sIt1 == rIt1->subRowsBegin() || sIt2 == rIt2->subRowsBegin())
                break;
            sIt1--;
            sIt2--;
        }

        if(verticalDir == 0)
        {
            ++rIt1;
            if(rIt1 == rbplace.coreRowsEnd())
                break;
        }
        else
        {
            if(rIt1 == rbplace.coreRowsBegin())
                break;
            else
                --rIt1;
        }
    }

    if(twoDim) {
        params.smplParams.algo = algo_save;
    }
}


bool divideWhitespace(const RBPSubRow &subrow, unsigned windowSize,
                      bool LR, SmallPlacementRow &smallRow, double currWSWidth, Point currWSLoc,
                      unsigned &movablesCounter, vector< vector<unsigned> > &movables,
                      vector<double> &whiteSpaceWidths, vector<Point> &whiteSpaceLocs,
                      unsigned wsPieceCap)
{
  vector<unsigned> cellIds(1, UINT_MAX);

  if(wsPieceCap == 0) return false;

  unsigned numSitesWS = static_cast<unsigned>(currWSWidth/smallRow.siteInterval);

  if(numSitesWS == 0) return false;

  unsigned wsPieces = std::min(numSitesWS, wsPieceCap);

//  make sure to include this entire piece of WS, not just
//  a fraction of it
//  wsPieces = std::min(wsPieces, windowSize - movablesCounter);

  unsigned sitesPerPiece = numSitesWS / wsPieces;
  unsigned sitesLeftOver = numSitesWS % wsPieces;

  for(unsigned i = 0; i < wsPieces; ++i)
  {
    unsigned sitesForThisPiece = (i >= wsPieces-sitesLeftOver) ? sitesPerPiece+1 : sitesPerPiece;

    movables.push_back(cellIds);
    whiteSpaceWidths.push_back(sitesForThisPiece*smallRow.siteInterval);
    whiteSpaceLocs.push_back(currWSLoc);
    if(LR)
    {
      currWSLoc.x += sitesForThisPiece*smallRow.siteInterval;
      smallRow.xMax += sitesForThisPiece*smallRow.siteInterval;
    }
    else
    {
      currWSLoc.x -= sitesForThisPiece*smallRow.siteInterval;
      smallRow.xMin -= sitesForThisPiece*smallRow.siteInterval;
    }
    ++movablesCounter;
    if(movablesCounter == windowSize)
    {
      return true;
    }
  }

  return false;
}
