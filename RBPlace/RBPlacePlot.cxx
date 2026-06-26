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

// by sadya
#ifdef _MSC_VER
#pragma warning(disable : 4786)
#endif

#include <cctype>
#include <cerrno>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>

#include <ABKCommon/abkassert.h>
#include <ABKCommon/infolines.h>
#include <ABKCommon/paramproc.h>
#include <ABKCommon/pathDelims.h>
#include <Ctainers/bitBoard.h>
#include <Geoms/bbox.h>
#include <HGraphWDims/hgWDims.h>
#include <PlaceEvals/pinPlEval.h>
#include <RBPlace/RBPlacePlot.h>
#include <RBPlace/RBPlacement.h>

#include <bitset>
#include <iostream>
#include <map>
#include <set>

using std::bitset;
using std::cin;
using std::cout;
using std::endl;
using std::map;
using std::ofstream;
using std::ostream;
using std::pair;
using std::set;
using std::setw;
using std::vector;
using std::string;

// TODO  Need to test histograph and sites map and row map etc.

// BEGIN RBPLACEPLOTTER CODE HERE
namespace Plotters {
RBPlacePlotter::Parameters::Parameters(void)
    : plotNodes(false),
      plotNodesNames(false),
      plotNets(false),
      plotSteinerNets(false),
      plotSiteMap(false),
      plotRows(false),
      plotFPOutlines(false),
      plotCutLines(false),
      plotFillerCells(true),
      oneColorPerGroup(false),
      xMin(-DBL_MAX),
      xMax(DBL_MAX),
      yMin(-DBL_MAX),
      yMax(DBL_MAX),
      FPOutlinesName(),
      CutLinesName() {}

RBPlacePlotter::Parameters::Parameters(const Parameters& orig)
    : plotNodes(orig.plotNodes),
      plotNodesNames(orig.plotNodesNames),
      plotNets(orig.plotNets),
      plotSteinerNets(orig.plotSteinerNets),
      plotSiteMap(orig.plotSiteMap),
      plotRows(orig.plotRows),
      plotFPOutlines(orig.plotFPOutlines),
      plotCutLines(orig.plotCutLines),
      plotFillerCells(orig.plotFillerCells),
      oneColorPerGroup(orig.oneColorPerGroup),
      xMin(orig.xMin),
      xMax(orig.xMax),
      yMin(orig.yMin),
      yMax(orig.yMax),
      FPOutlinesName(orig.FPOutlinesName),
      CutLinesName(orig.CutLinesName) {}

RBPlacePlotter::Parameters::Parameters(int argc, const char** argv,
                                       bool switchBoxEnabled)
    : plotNodes(false),
      plotNodesNames(false),
      plotNets(false),
      plotSteinerNets(false),
      plotSiteMap(false),
      plotRows(false),
      plotFPOutlines(false),
      plotCutLines(false),
      plotFillerCells(true),
      oneColorPerGroup(false),
      xMin(-DBL_MAX),
      xMax(DBL_MAX),
      yMin(-DBL_MAX),
      yMax(DBL_MAX) {
  if (switchBoxEnabled) {
    BoolParam plot_Nodes("dpPlotNodes", argc, argv);
    plotNodes = plotNodes || plot_Nodes.found();

    BoolParam plot_Nodes_Names("dpPlotNodesNames", argc, argv);
    plotNodesNames = plotNodesNames || plot_Nodes_Names.found();

    BoolParam plot_Nets("dpPlotNets", argc, argv);
    plotNets = plotNets || plot_Nets.found();

    BoolParam plot_Site_Map("dpPlotSiteMap", argc, argv);
    plotSiteMap = plotSiteMap || plot_Site_Map.found();

    BoolParam plot_Rows("dpPlotRows", argc, argv);
    plotRows = plotRows || plot_Rows.found();
  }

  DoubleParam xmin("xmin", argc, argv);
  if (xmin.found()) xMin = xmin;

  DoubleParam xmax("xmax", argc, argv);
  if (xmax.found()) xMax = xmax;

  DoubleParam ymin("ymin", argc, argv);
  if (ymin.found()) yMin = ymin;

  DoubleParam ymax("ymax", argc, argv);
  if (ymax.found()) yMax = ymax;
}

RBPlacePlotter::RBPlacePlotter(RBPlacePlotter::Parameters modes,
                               std::string baseFileName,
                               const RBPlacement& rbp)
    : _modes(modes),
      _baseFileName(baseFileName),
      _rbp(rbp),
      _plotBox(rbp.getBBox()),
      _file_ct(0) {}

std::ostream& operator<<(std::ostream& os, const RBPlacePlotter& rbp) {
  rbp.writeGnuplotFileHeader(os);
  rbp.writeGnuplotDataFilesNameLabels(os);
  rbp.writeGnuplotPlotCommand(os);
  rbp.writeGnuplotDataFilesContents(os);
  rbp.writeGnuplotFileFooter(os);
  return os;
}

void RBPlacePlotter::update(void) {
  setFilesAndColors();
  return;
}

void RBPlacePlotter::setFilesAndColors(void) {
  _file_ct = getFileCount();

  // if(! user defined colors) //TODO
  assignStandardColors();  // TODO
  return;
}

void RBPlacePlotter::writeGnuplotDataFilesNameLabels(std::ostream& os) const {
  if (_modes.plotNodesNames) {
    for (vector<nodeSet>::const_iterator it = _nodes2Label.begin();
         it != _nodes2Label.end(); ++it)
      writeGnuplotNodeNameLabels(os, *it);
  }
}

void RBPlacePlotter::writeGnuplotDataFilesContents(std::ostream& os) const {
  if (_modes.plotNodes) {
    for (vector<nodeSet>::const_iterator it = _nodeGroups.begin();
         it != _nodeGroups.end(); ++it) {
      writeGnuplotSmallMovableNodes(os, *it);
      writeGnuplotSmallFixedNodes(os, *it);
      writeGnuplotLargeMovableNodes(os, *it);
      writeGnuplotLargeFixedNodes(os, *it);
    }
  }

  if (_modes.plotNets) {
    if (_modes.plotSteinerNets) {
      for (vector<nodeSet>::const_iterator it = _nodes4Nets.begin();
           it != _nodes4Nets.end(); ++it) {
        writeGnuplotSteinerNets(os, *it);
        os << "EOF" << endl << endl;
      }
    } else {
      for (vector<nodeSet>::const_iterator it = _nodes4Nets.begin();
           it != _nodes4Nets.end(); ++it) {
        writeGnuplotNets(os, *it);
        os << "EOF" << endl << endl;
      }
    }
  }

  if (_modes.plotSiteMap) writeGnuplotSiteMap(os);

  if (_modes.plotRows) writeGnuplotRows(os);
}

void RBPlacePlotter::writeGnuplotFileHeader(std::ostream& os) const {
  os << "#Use this file as a script for gnuplot" << endl;
  os << "#(See http://www.gnuplot.info/ for details)" << endl;
  os << TimeStamp() << User() << Platform() << endl;
  unsigned oldprecision = os.precision(4);
  os << "set title \" " << _baseFileName << " HPWL= " << _rbp.evalHPWL(true)
     << ", #Cells= " << _rbp.getNumCells()
     << ", #Nets= " << (_rbp.getHGraph().getNumEdges())
     << " \" font \"Times , 22\"" << endl
     << endl;
  os.precision(oldprecision);
  os << "set nokey" << endl;

  os << "#   Uncomment these two lines starting with \"set\"" << endl;
  os << "#   to save an EPS file for inclusion into a latex document" << endl;
  os << "# set terminal postscript eps color solid 20" << endl;
  os << "# set output \"" << _baseFileName << ".eps\"" << endl << endl;

  os << "#   Uncomment these two lines starting with \"set\"" << endl;
  os << "#   to save a PS file for printing" << endl;
  os << "# set terminal postscript portrait color solid 20" << endl;
  os << "# set output \"" << _baseFileName << ".ps\"" << endl << endl;
}

void RBPlacePlotter::writeGnuplotFileFooter(std::ostream& os) const {
  os << "pause -1 'Press any key' " << endl;
}

void RBPlacePlotter::writeGnuplotNodeNameLabels(std::ostream& os,
                                                const nodeSet& nodes) const {
  for (vector<nodeSet>::const_iterator it = _nodeGroups.begin();
       it != _nodeGroups.end(); ++it) {
    const nodeSet& nodes = *it;
    for (size_t i = 0; i < nodes.size(); i++) {
      unsigned index = nodes[i];
      os << "set label \"" << _rbp.getHGraph().getNodeNameByIndex(index)
         << "\" at "
         << (((_rbp.getPlacement())[index]).x) +
                _rbp.getHGraph().getNodeWidth(index) / 4
         << " , "
         << _rbp.getPlacement()[index].y +
                _rbp.getHGraph().getNodeHeight(index) / 3
         << endl;
    }
  }
  os << endl << endl;
}

void RBPlacePlotter::writeGnuplotPlotCommand(std::ostream& os) const {
  os << " plot[";
  if (_modes.xMin != -DBL_MAX) os << _modes.xMin;
  os << ":";
  if (_modes.xMax != DBL_MAX) os << _modes.xMax;

  os << "][";

  if (_modes.yMin != -DBL_MAX) os << _modes.yMin;
  os << ":";
  if (_modes.yMax != DBL_MAX) os << _modes.yMax;

  os << "]";

  //        if(_colors.size() != _file_ct )
  //        {
  //            cout<<"There are "<<_colors.size()<<" colors and "<<_file_ct<<"
  //            files" <<endl; abkfatal(_colors.size() != _file_ct, "Colors and
  //            Files do not match, Internal Capo Error, Please report this
  //            bug");
  //        }

  if (_colors.size() > 0) {
    for (size_t index = 0; index < _colors.size() - 1; ++index) {
      os << " '-' w l lt " << _colors[index] << ",";
    }
    os << " '-' w l lt " << _colors.back();
  }
  if (_modes.plotCutLines) {
    if (_file_ct > 0) os << ",";
    os << " '" << _modes.CutLinesName << "' w l lt " << Plotters::red
       << " lw 2";
  }
  if (_modes.plotFPOutlines) {
    if (_file_ct > 0 || _modes.plotCutLines) os << ",";
    os << " '" << _modes.FPOutlinesName << "' w l lt " << Plotters::red;
  }
  os << endl << endl;
}

void RBPlacePlotter::writeGnuplotNets(std::ostream& os,
                                      const nodeSet& nodes) const {
  // the hyperedges info goes inside
  vector<bool> seenNets(_rbp.getHGraph().getNumEdges(), false);
  std::string cellName;
  unsigned netDegree = 0;
  Point lastRowPoint;
  double xCentre = 0, yCentre = 0;
  unsigned i = 0;
  double rowHeight = 0;
  double nodeWidth = 0;
  double rowLength = 1;
  vector<double> xcoord;
  vector<double> ycoord;

  Point minRange;
  Point maxRange;

  minRange.x = _modes.xMin;
  maxRange.x = _modes.xMax;
  minRange.y = _modes.yMin;
  maxRange.y = _modes.yMax;

  for (unsigned k = 0; k < nodes.size(); ++k) {
    unsigned index = nodes[k];
    const HGFNode& node = _rbp.getHGraph().getNodeByIdx(index);

    for (itHGFEdgeLocal e = node.edgesBegin(); e != node.edgesEnd(); e++) {
      const HGFEdge& edge = (**e);
      unsigned edgeId = edge.getIndex();
      if (!seenNets[edgeId]) {
        seenNets[edgeId] = true;
        netDegree = edge.getDegree();
        os << "#NetDegree : " << netDegree << endl;

        xcoord.resize(netDegree);
        ycoord.resize(netDegree);
        unsigned nodeOffset = 0;

        itHGFNodeLocal v;
        unsigned j = 0;

        unsigned withinRange = 0;
        unsigned status = 0;
        xCentre = 0;
        yCentre = 0;
        double ycordRow = 0;

        // for all nodes
        for (v = edge.nodesBegin(), j = 0; v != edge.nodesEnd();
             v++, nodeOffset++) {
          cellName = _rbp.getHGraph().getNodeNameByIndex((*v)->getIndex());
          os << "#" << setw(10) << cellName << " O " << endl;
          //                    HGFNode& node =
          //                    _rbp.getHGraph().getNodeByName(cellName.c_str());
          //                    i=node.getIndex();
          i = (*v)->getIndex();
          Point tempRowPoint(_rbp.getPlacement()[i].x,
                             _rbp.getPlacement()[i].y);
          lastRowPoint = tempRowPoint;

          if (j == 0)  // ycordRow gives the ycord of first cell of row
            ycordRow = _rbp.getPlacement()[i]
                           .y;  // used to determine if cells in 1 row

          Point pOffset = _rbp.getHGraph().nodesOnEdgePinOffset(
              nodeOffset, edgeId, _rbp.getPlacement().getOrient(i));

          ycoord[j] = _rbp.getPlacement()[i].y + pOffset.y;
          xcoord[j] = _rbp.getPlacement()[i].x + pOffset.x;
          rowHeight = _rbp.getHGraph().getNodeHeight(i);
          nodeWidth = _rbp.getHGraph().getNodeWidth(i);

          if ((_rbp.getPlacement()[i].x > minRange.x &&
               _rbp.getPlacement()[i].x < maxRange.x) &&
              (_rbp.getPlacement()[i].y > minRange.y &&
               _rbp.getPlacement()[i].y < maxRange.y))
            withinRange = 1;  // atleast 1 cell of Hedge is within range
          // so plot the whole net
          if (ycordRow != _rbp.getPlacement()[i].y) {
            status = 1;  // all cells of Hedge are not in single row
          }

          xCentre += xcoord[j];
          yCentre += ycoord[j];
          j++;
        }
        if (withinRange == 1) {
          if (status == 1) {
            xCentre /= netDegree;
            yCentre /= netDegree;
          } else {
            const RBPCoreRow* rowName;
            rowName = _rbp.findCoreRow(lastRowPoint);
            if (rowName != NULL)
              rowLength = fabs(rowName->getXEnd() - rowName->getXStart());
            double maxXcoord = 0;
            double minXcoord = 0;

            minXcoord = *std::min_element(xcoord.begin(), xcoord.end());
            maxXcoord = *std::max_element(xcoord.begin(), xcoord.end());
            double ydifference = maxXcoord - minXcoord;
            double rd = RandomDouble(-1, 1);

            int starOrient = 0;
            starOrient = ((rd > 0) ? 1 : -1);  // get orientation of star point

            rd = RandomDouble(-1, 1);
            ydifference =
                (1 - exp(0 - ydifference * ydifference / (rowLength))) *
                rowHeight * 0.4 * starOrient;
            yCentre = ycoord[0] + ydifference;
            xCentre = (maxXcoord + minXcoord) / 2 + rd * nodeWidth * 0.25;
          }
          if (status == 0 || netDegree != 2) {
            if (fabs(xCentre) < 1e+306 && fabs(yCentre) < 1e+306) {
              os << setw(10) << xCentre << "  " << setw(10)
                 << yCentre + 0.025 * rowHeight << "  " << endl;
              os << setw(10) << xCentre + 0.025 * rowHeight << "  " << setw(10)
                 << yCentre << "  " << endl;
              os << setw(10) << xCentre << "  " << setw(10)
                 << yCentre - 0.025 * rowHeight << "  " << endl;
              os << setw(10) << xCentre - 0.025 * rowHeight << "  " << setw(10)
                 << yCentre << "  " << endl;
              os << setw(10) << xCentre << "  " << setw(10)
                 << yCentre + 0.025 * rowHeight << "  " << endl
                 << endl;
            }
          }

          for (i = 0; i < netDegree; ++i) {
            if (fabs(xcoord[i]) < DBL_MAX && fabs(ycoord[i]) < DBL_MAX &&
                fabs(xCentre) < DBL_MAX && fabs(yCentre) < DBL_MAX) {
              os << setw(10) << xCentre << "  " << setw(10) << yCentre << "  "
                 << endl;
              os << setw(10) << xcoord[i] << "  " << setw(10) << ycoord[i]
                 << "  " << endl;
              os << setw(10) << xCentre << "  " << setw(10) << yCentre << "  "
                 << endl
                 << endl;
            }
          }
        }
      }
    }
  }
}

void RBPlacePlotter::writeGnuplotSteinerNets(std::ostream& os,
                                             const nodeSet& nodes) const {
  // the hyperedges info goes inside
  BitBoard seenNets(_rbp.getHGraph().getNumEdges());

  std::vector<pair<Point, Point> > edges;
  bool usePinLocs = true;
  bool cleanup = true;

  for (unsigned k = 0; k < nodes.size(); ++k) {
    unsigned index = nodes[k];
    const HGFNode& node = _rbp.getHGraph().getNodeByIdx(index);

    for (itHGFEdgeLocal e = node.edgesBegin(); e != node.edgesEnd(); ++e) {
      unsigned edgeId = (*e)->getIndex();
      seenNets.setBit(edgeId);
    }
  }

  const auto& netsToPlot = seenNets.getIndicesOfSetBits();
  unsigned cutoff1 = netsToPlot.size() / 3;
  unsigned cutoff2 = (2 * netsToPlot.size()) / 3;

  for (unsigned i = 0; i < netsToPlot.size(); ++i) {
    if (i == cutoff1 || i == cutoff2) {
      os << endl << "EOF" << endl;
    }

    pinPlEval::getBestSteinerTree(_rbp.getPlacement(), _rbp.getHGraph(),
                                  netsToPlot[i], usePinLocs, cleanup, edges);

    for (unsigned j = 0; j < edges.size(); ++j) {
      os << edges[j].first.x << " " << edges[j].first.y << endl
         << edges[j].second.x << " " << edges[j].second.y << endl
         << endl;
    }
  }
}

void RBPlacePlotter::writeGnuplotLargeFixedNodes(std::ostream& os,
                                                 const nodeSet& nodes) const {
  double cellWidth = 0;
  double cellHeight = 0;
  unsigned j;
  unsigned i = 0;
  bool haveAny = false;
  for (j = 0; j < nodes.size(); j++) {
    i = nodes[j];
    const Orient& nodeOrient = _rbp.getPlacement().getOrient(i);
    unsigned angle = nodeOrient.getAngle();
    bool notRotated = angle == 0 || angle == 180;
    cellWidth = notRotated ? _rbp.getHGraph().getNodeWidth(i)
                           : _rbp.getHGraph().getNodeHeight(i);
    cellHeight = notRotated ? _rbp.getHGraph().getNodeHeight(i)
                            : _rbp.getHGraph().getNodeWidth(i);

    double maxCRHeight = _rbp.getMaxHeightCoreRow();
    if (((cellHeight - maxCRHeight) > 1e-6) && _rbp.isFixed(i)) {
      haveAny = true;
      double offset = 0;
      offset =
          ((cellWidth > cellHeight) ? 0.15 * cellHeight : 0.15 * cellWidth);
      os << "# Cell Name " << _rbp.getHGraph().getNodeNameByIndex(i) << endl;

      if (fabs(_rbp.getPlacement()[i].x) < 1e+306 &&
          fabs(_rbp.getPlacement()[i].y) < 1e+306) {
        if (nodeOrient == Orient("S")) {
          os << setw(10) << _rbp.getPlacement()[i].x << "  " << setw(10)
             << _rbp.getPlacement()[i].y + offset << endl;
          os << setw(10) << _rbp.getPlacement()[i].x + 3. * offset << "  "
             << setw(10) << _rbp.getPlacement()[i].y << endl;
        } else if (nodeOrient == Orient("FE")) {
          os << setw(10) << _rbp.getPlacement()[i].x << "  " << setw(10)
             << _rbp.getPlacement()[i].y + 3. * offset << endl;
          os << setw(10) << _rbp.getPlacement()[i].x + offset << "  "
             << setw(10) << _rbp.getPlacement()[i].y << endl;
        } else {
          os << setw(10) << _rbp.getPlacement()[i].x << "  " << setw(10)
             << _rbp.getPlacement()[i].y << endl;
        }

        if (nodeOrient == Orient("E")) {
          os << setw(10) << _rbp.getPlacement()[i].x + cellWidth - offset
             << "  " << setw(10) << _rbp.getPlacement()[i].y << endl;
          os << setw(10) << _rbp.getPlacement()[i].x + cellWidth << "  "
             << setw(10) << _rbp.getPlacement()[i].y + 3. * offset << endl;
        } else if (nodeOrient == Orient("FS")) {
          os << setw(10) << _rbp.getPlacement()[i].x + cellWidth - 3. * offset
             << "  " << setw(10) << _rbp.getPlacement()[i].y << endl;
          os << setw(10) << _rbp.getPlacement()[i].x + cellWidth << "  "
             << setw(10) << _rbp.getPlacement()[i].y + offset << endl;
        } else {
          os << setw(10) << _rbp.getPlacement()[i].x + cellWidth << "  "
             << setw(10) << _rbp.getPlacement()[i].y << endl;
        }

        if (nodeOrient == Orient("N")) {
          os << setw(10) << _rbp.getPlacement()[i].x + cellWidth << "  "
             << setw(10) << _rbp.getPlacement()[i].y + cellHeight - offset
             << endl;
          os << setw(10) << _rbp.getPlacement()[i].x + cellWidth - 3. * offset
             << "  " << setw(10) << _rbp.getPlacement()[i].y + cellHeight
             << endl;
        } else if (nodeOrient == Orient("FW")) {
          os << setw(10) << _rbp.getPlacement()[i].x + cellWidth << "  "
             << setw(10) << _rbp.getPlacement()[i].y + cellHeight - 3. * offset
             << endl;
          os << setw(10) << _rbp.getPlacement()[i].x + cellWidth - offset
             << "  " << setw(10) << _rbp.getPlacement()[i].y + cellHeight
             << endl;
        } else {
          os << setw(10) << _rbp.getPlacement()[i].x + cellWidth << "  "
             << setw(10) << _rbp.getPlacement()[i].y + cellHeight << endl;
        }

        if (nodeOrient == Orient("W")) {
          os << setw(10) << _rbp.getPlacement()[i].x + offset << "  "
             << setw(10) << _rbp.getPlacement()[i].y + cellHeight << endl;
          os << setw(10) << _rbp.getPlacement()[i].x << "  " << setw(10)
             << _rbp.getPlacement()[i].y + cellHeight - 3. * offset << endl;
        } else if (nodeOrient == Orient("FN")) {
          os << setw(10) << _rbp.getPlacement()[i].x + 3. * offset << "  "
             << setw(10) << _rbp.getPlacement()[i].y + cellHeight << endl;
          os << setw(10) << _rbp.getPlacement()[i].x << "  " << setw(10)
             << _rbp.getPlacement()[i].y + cellHeight - offset << endl;
        } else {
          os << setw(10) << _rbp.getPlacement()[i].x << "  " << setw(10)
             << _rbp.getPlacement()[i].y + cellHeight << endl;
        }

        if (nodeOrient == Orient("S")) {
          os << setw(10) << _rbp.getPlacement()[i].x << "  " << setw(10)
             << _rbp.getPlacement()[i].y + offset << endl;
        } else if (nodeOrient == Orient("FE")) {
          os << setw(10) << _rbp.getPlacement()[i].x << "  " << setw(10)
             << _rbp.getPlacement()[i].y + 3. * offset << endl;
        } else {
          os << setw(10) << _rbp.getPlacement()[i].x << "  " << setw(10)
             << _rbp.getPlacement()[i].y << endl;
        }

        os << endl << endl;

        // do the inner cell
        double innerOffset = 0.667 * offset;
        offset *= 0.8;

        if (nodeOrient == Orient("S")) {
          os << setw(10) << _rbp.getPlacement()[i].x + innerOffset << "  "
             << setw(10) << _rbp.getPlacement()[i].y + innerOffset + offset
             << endl;
          os << setw(10) << _rbp.getPlacement()[i].x + innerOffset + 3. * offset
             << "  " << setw(10) << _rbp.getPlacement()[i].y + innerOffset
             << endl;
        } else if (nodeOrient == Orient("FE")) {
          os << setw(10) << _rbp.getPlacement()[i].x + innerOffset << "  "
             << setw(10) << _rbp.getPlacement()[i].y + innerOffset + 3. * offset
             << endl;
          os << setw(10) << _rbp.getPlacement()[i].x + innerOffset + offset
             << "  " << setw(10) << _rbp.getPlacement()[i].y + innerOffset
             << endl;
        } else {
          os << setw(10) << _rbp.getPlacement()[i].x + innerOffset << "  "
             << setw(10) << _rbp.getPlacement()[i].y + innerOffset << endl;
        }

        if (nodeOrient == Orient("E")) {
          os << setw(10)
             << _rbp.getPlacement()[i].x - innerOffset + cellWidth - offset
             << "  " << setw(10) << _rbp.getPlacement()[i].y + innerOffset
             << endl;
          os << setw(10) << _rbp.getPlacement()[i].x - innerOffset + cellWidth
             << "  " << setw(10)
             << _rbp.getPlacement()[i].y + innerOffset + 3. * offset << endl;
        } else if (nodeOrient == Orient("FS")) {
          os << setw(10)
             << _rbp.getPlacement()[i].x - innerOffset + cellWidth - 3. * offset
             << "  " << setw(10) << _rbp.getPlacement()[i].y + innerOffset
             << endl;
          os << setw(10) << _rbp.getPlacement()[i].x - innerOffset + cellWidth
             << "  " << setw(10)
             << _rbp.getPlacement()[i].y + innerOffset + offset << endl;
        } else {
          os << setw(10) << _rbp.getPlacement()[i].x - innerOffset + cellWidth
             << "  " << setw(10) << _rbp.getPlacement()[i].y + innerOffset
             << endl;
        }

        if (nodeOrient == Orient("N")) {
          os << setw(10) << _rbp.getPlacement()[i].x - innerOffset + cellWidth
             << "  " << setw(10)
             << _rbp.getPlacement()[i].y - innerOffset + cellHeight - offset
             << endl;
          os << setw(10)
             << _rbp.getPlacement()[i].x - innerOffset + cellWidth - 3. * offset
             << "  " << setw(10)
             << _rbp.getPlacement()[i].y - innerOffset + cellHeight << endl;
        } else if (nodeOrient == Orient("FW")) {
          os << setw(10) << _rbp.getPlacement()[i].x - innerOffset + cellWidth
             << "  " << setw(10)
             << _rbp.getPlacement()[i].y - innerOffset + cellHeight -
                    3. * offset
             << endl;
          os << setw(10)
             << _rbp.getPlacement()[i].x - innerOffset + cellWidth - offset
             << "  " << setw(10)
             << _rbp.getPlacement()[i].y - innerOffset + cellHeight << endl;
        } else {
          os << setw(10) << _rbp.getPlacement()[i].x - innerOffset + cellWidth
             << "  " << setw(10)
             << _rbp.getPlacement()[i].y - innerOffset + cellHeight << endl;
        }

        if (nodeOrient == Orient("W")) {
          os << setw(10) << _rbp.getPlacement()[i].x + innerOffset + offset
             << "  " << setw(10)
             << _rbp.getPlacement()[i].y - innerOffset + cellHeight << endl;
          os << setw(10) << _rbp.getPlacement()[i].x + innerOffset << "  "
             << setw(10)
             << _rbp.getPlacement()[i].y - innerOffset + cellHeight -
                    3. * offset
             << endl;
        } else if (nodeOrient == Orient("FN")) {
          os << setw(10) << _rbp.getPlacement()[i].x + innerOffset + 3. * offset
             << "  " << setw(10)
             << _rbp.getPlacement()[i].y - innerOffset + cellHeight << endl;
          os << setw(10) << _rbp.getPlacement()[i].x + innerOffset << "  "
             << setw(10)
             << _rbp.getPlacement()[i].y - innerOffset + cellHeight - offset
             << endl;
        } else {
          os << setw(10) << _rbp.getPlacement()[i].x + innerOffset << "  "
             << setw(10) << _rbp.getPlacement()[i].y - innerOffset + cellHeight
             << endl;
        }

        if (nodeOrient == Orient("S")) {
          os << setw(10) << _rbp.getPlacement()[i].x + innerOffset << "  "
             << setw(10) << _rbp.getPlacement()[i].y + innerOffset + offset
             << endl;
        } else if (nodeOrient == Orient("FE")) {
          os << setw(10) << _rbp.getPlacement()[i].x + innerOffset << "  "
             << setw(10) << _rbp.getPlacement()[i].y + innerOffset + 3. * offset
             << endl;
        } else {
          os << setw(10) << _rbp.getPlacement()[i].x + innerOffset << "  "
             << setw(10) << _rbp.getPlacement()[i].y + innerOffset << endl;
        }

        os << endl << endl;
      }
    }
  }

  if (!haveAny) os << "0" << endl << endl;

  os << "EOF" << endl << endl;
}

void RBPlacePlotter::writeGnuplotLargeMovableNodes(std::ostream& os,
                                                   const nodeSet& nodes) const {
  double cellWidth = 0;
  double cellHeight = 0;
  unsigned j;
  unsigned i = 0;
  bool haveAny = false;
  for (j = 0; j < nodes.size(); j++) {
    i = nodes[j];
    const Orient& nodeOrient = _rbp.getPlacement().getOrient(i);
    unsigned angle = nodeOrient.getAngle();
    bool notRotated = angle == 0 || angle == 180;
    cellWidth = notRotated ? _rbp.getHGraph().getNodeWidth(i)
                           : _rbp.getHGraph().getNodeHeight(i);
    cellHeight = notRotated ? _rbp.getHGraph().getNodeHeight(i)
                            : _rbp.getHGraph().getNodeWidth(i);

    double maxCRHeight = _rbp.getMaxHeightCoreRow();
    if (((cellHeight - maxCRHeight) > 1e-6) && !_rbp.isFixed(i)) {
      haveAny = true;
      double offset = 0;
      offset =
          ((cellWidth > cellHeight) ? 0.15 * cellHeight : 0.15 * cellWidth);
      os << "# Cell Name " << _rbp.getHGraph().getNodeNameByIndex(i) << endl;

      if (fabs(_rbp.getPlacement()[i].x) < 1e+306 &&
          fabs(_rbp.getPlacement()[i].y) < 1e+306) {
        if (nodeOrient == Orient("S")) {
          os << setw(10) << _rbp.getPlacement()[i].x << "  " << setw(10)
             << _rbp.getPlacement()[i].y + offset << endl;
          os << setw(10) << _rbp.getPlacement()[i].x + 3. * offset << "  "
             << setw(10) << _rbp.getPlacement()[i].y << endl;
        } else if (nodeOrient == Orient("FE")) {
          os << setw(10) << _rbp.getPlacement()[i].x << "  " << setw(10)
             << _rbp.getPlacement()[i].y + 3. * offset << endl;
          os << setw(10) << _rbp.getPlacement()[i].x + offset << "  "
             << setw(10) << _rbp.getPlacement()[i].y << endl;
        } else {
          os << setw(10) << _rbp.getPlacement()[i].x << "  " << setw(10)
             << _rbp.getPlacement()[i].y << endl;
        }

        if (nodeOrient == Orient("E")) {
          os << setw(10) << _rbp.getPlacement()[i].x + cellWidth - offset
             << "  " << setw(10) << _rbp.getPlacement()[i].y << endl;
          os << setw(10) << _rbp.getPlacement()[i].x + cellWidth << "  "
             << setw(10) << _rbp.getPlacement()[i].y + 3. * offset << endl;
        } else if (nodeOrient == Orient("FS")) {
          os << setw(10) << _rbp.getPlacement()[i].x + cellWidth - 3. * offset
             << "  " << setw(10) << _rbp.getPlacement()[i].y << endl;
          os << setw(10) << _rbp.getPlacement()[i].x + cellWidth << "  "
             << setw(10) << _rbp.getPlacement()[i].y + offset << endl;
        } else {
          os << setw(10) << _rbp.getPlacement()[i].x + cellWidth << "  "
             << setw(10) << _rbp.getPlacement()[i].y << endl;
        }

        if (nodeOrient == Orient("N")) {
          os << setw(10) << _rbp.getPlacement()[i].x + cellWidth << "  "
             << setw(10) << _rbp.getPlacement()[i].y + cellHeight - offset
             << endl;
          os << setw(10) << _rbp.getPlacement()[i].x + cellWidth - 3. * offset
             << "  " << setw(10) << _rbp.getPlacement()[i].y + cellHeight
             << endl;
        } else if (nodeOrient == Orient("FW")) {
          os << setw(10) << _rbp.getPlacement()[i].x + cellWidth << "  "
             << setw(10) << _rbp.getPlacement()[i].y + cellHeight - 3. * offset
             << endl;
          os << setw(10) << _rbp.getPlacement()[i].x + cellWidth - offset
             << "  " << setw(10) << _rbp.getPlacement()[i].y + cellHeight
             << endl;
        } else {
          os << setw(10) << _rbp.getPlacement()[i].x + cellWidth << "  "
             << setw(10) << _rbp.getPlacement()[i].y + cellHeight << endl;
        }

        if (nodeOrient == Orient("W")) {
          os << setw(10) << _rbp.getPlacement()[i].x + offset << "  "
             << setw(10) << _rbp.getPlacement()[i].y + cellHeight << endl;
          os << setw(10) << _rbp.getPlacement()[i].x << "  " << setw(10)
             << _rbp.getPlacement()[i].y + cellHeight - 3. * offset << endl;
        } else if (nodeOrient == Orient("FN")) {
          os << setw(10) << _rbp.getPlacement()[i].x + 3. * offset << "  "
             << setw(10) << _rbp.getPlacement()[i].y + cellHeight << endl;
          os << setw(10) << _rbp.getPlacement()[i].x << "  " << setw(10)
             << _rbp.getPlacement()[i].y + cellHeight - offset << endl;
        } else {
          os << setw(10) << _rbp.getPlacement()[i].x << "  " << setw(10)
             << _rbp.getPlacement()[i].y + cellHeight << endl;
        }

        if (nodeOrient == Orient("S")) {
          os << setw(10) << _rbp.getPlacement()[i].x << "  " << setw(10)
             << _rbp.getPlacement()[i].y + offset << endl;
        } else if (nodeOrient == Orient("FE")) {
          os << setw(10) << _rbp.getPlacement()[i].x << "  " << setw(10)
             << _rbp.getPlacement()[i].y + 3. * offset << endl;
        } else {
          os << setw(10) << _rbp.getPlacement()[i].x << "  " << setw(10)
             << _rbp.getPlacement()[i].y << endl;
        }

        os << endl << endl;

        // do the inner cell
        double innerOffset = 0.667 * offset;
        offset *= 0.8;

        if (nodeOrient == Orient("S")) {
          os << setw(10) << _rbp.getPlacement()[i].x + innerOffset << "  "
             << setw(10) << _rbp.getPlacement()[i].y + innerOffset + offset
             << endl;
          os << setw(10) << _rbp.getPlacement()[i].x + innerOffset + 3. * offset
             << "  " << setw(10) << _rbp.getPlacement()[i].y + innerOffset
             << endl;
        } else if (nodeOrient == Orient("FE")) {
          os << setw(10) << _rbp.getPlacement()[i].x + innerOffset << "  "
             << setw(10) << _rbp.getPlacement()[i].y + innerOffset + 3. * offset
             << endl;
          os << setw(10) << _rbp.getPlacement()[i].x + innerOffset + offset
             << "  " << setw(10) << _rbp.getPlacement()[i].y + innerOffset
             << endl;
        } else {
          os << setw(10) << _rbp.getPlacement()[i].x + innerOffset << "  "
             << setw(10) << _rbp.getPlacement()[i].y + innerOffset << endl;
        }

        if (nodeOrient == Orient("E")) {
          os << setw(10)
             << _rbp.getPlacement()[i].x - innerOffset + cellWidth - offset
             << "  " << setw(10) << _rbp.getPlacement()[i].y + innerOffset
             << endl;
          os << setw(10) << _rbp.getPlacement()[i].x - innerOffset + cellWidth
             << "  " << setw(10)
             << _rbp.getPlacement()[i].y + innerOffset + 3. * offset << endl;
        } else if (nodeOrient == Orient("FS")) {
          os << setw(10)
             << _rbp.getPlacement()[i].x - innerOffset + cellWidth - 3. * offset
             << "  " << setw(10) << _rbp.getPlacement()[i].y + innerOffset
             << endl;
          os << setw(10) << _rbp.getPlacement()[i].x - innerOffset + cellWidth
             << "  " << setw(10)
             << _rbp.getPlacement()[i].y + innerOffset + offset << endl;
        } else {
          os << setw(10) << _rbp.getPlacement()[i].x - innerOffset + cellWidth
             << "  " << setw(10) << _rbp.getPlacement()[i].y + innerOffset
             << endl;
        }

        if (nodeOrient == Orient("N")) {
          os << setw(10) << _rbp.getPlacement()[i].x - innerOffset + cellWidth
             << "  " << setw(10)
             << _rbp.getPlacement()[i].y - innerOffset + cellHeight - offset
             << endl;
          os << setw(10)
             << _rbp.getPlacement()[i].x - innerOffset + cellWidth - 3. * offset
             << "  " << setw(10)
             << _rbp.getPlacement()[i].y - innerOffset + cellHeight << endl;
        } else if (nodeOrient == Orient("FW")) {
          os << setw(10) << _rbp.getPlacement()[i].x - innerOffset + cellWidth
             << "  " << setw(10)
             << _rbp.getPlacement()[i].y - innerOffset + cellHeight -
                    3. * offset
             << endl;
          os << setw(10)
             << _rbp.getPlacement()[i].x - innerOffset + cellWidth - offset
             << "  " << setw(10)
             << _rbp.getPlacement()[i].y - innerOffset + cellHeight << endl;
        } else {
          os << setw(10) << _rbp.getPlacement()[i].x - innerOffset + cellWidth
             << "  " << setw(10)
             << _rbp.getPlacement()[i].y - innerOffset + cellHeight << endl;
        }

        if (nodeOrient == Orient("W")) {
          os << setw(10) << _rbp.getPlacement()[i].x + innerOffset + offset
             << "  " << setw(10)
             << _rbp.getPlacement()[i].y - innerOffset + cellHeight << endl;
          os << setw(10) << _rbp.getPlacement()[i].x + innerOffset << "  "
             << setw(10)
             << _rbp.getPlacement()[i].y - innerOffset + cellHeight -
                    3. * offset
             << endl;
        } else if (nodeOrient == Orient("FN")) {
          os << setw(10) << _rbp.getPlacement()[i].x + innerOffset + 3. * offset
             << "  " << setw(10)
             << _rbp.getPlacement()[i].y - innerOffset + cellHeight << endl;
          os << setw(10) << _rbp.getPlacement()[i].x + innerOffset << "  "
             << setw(10)
             << _rbp.getPlacement()[i].y - innerOffset + cellHeight - offset
             << endl;
        } else {
          os << setw(10) << _rbp.getPlacement()[i].x + innerOffset << "  "
             << setw(10) << _rbp.getPlacement()[i].y - innerOffset + cellHeight
             << endl;
        }

        if (nodeOrient == Orient("S")) {
          os << setw(10) << _rbp.getPlacement()[i].x + innerOffset << "  "
             << setw(10) << _rbp.getPlacement()[i].y + innerOffset + offset
             << endl;
        } else if (nodeOrient == Orient("FE")) {
          os << setw(10) << _rbp.getPlacement()[i].x + innerOffset << "  "
             << setw(10) << _rbp.getPlacement()[i].y + innerOffset + 3. * offset
             << endl;
        } else {
          os << setw(10) << _rbp.getPlacement()[i].x + innerOffset << "  "
             << setw(10) << _rbp.getPlacement()[i].y + innerOffset << endl;
        }

        os << endl << endl;
      }
    }
  }

  if (!haveAny) os << "0" << endl << endl;

  os << "EOF" << endl << endl;
}

void RBPlacePlotter::writeGnuplotSmallFixedNodes(std::ostream& os,
                                                 const nodeSet& nodes) const {
  bool haveAny = false;
  vector<vector<bool> >*horizPixelBuf = 0, *vertPixelBuf = 0;
  double cellWidth = 0., cellHeight = 0.;
  unsigned i = 0, j;
  double unitsperx = 0.0, unitspery = 0.0;
  double minx = DBL_MAX, miny = DBL_MAX, maxx = -DBL_MAX, maxy = -DBL_MAX;
  if (_rbp.getParameters().adaptivePlotting) {
    for (j = 0; j < nodes.size(); j++) {
      i = nodes[j];
      if (fabs(_rbp.getPlacement()[i].x) >= 1e+306 ||
          fabs(_rbp.getPlacement()[i].y) >= 1e+306)
        continue;
      const Orient& nodeOrient = _rbp.getPlacement().getOrient(i);
      unsigned angle = nodeOrient.getAngle();
      bool notRotated = angle == 0 || angle == 180;
      cellWidth = notRotated ? _rbp.getHGraph().getNodeWidth(i)
                             : _rbp.getHGraph().getNodeHeight(i);
      cellHeight = notRotated ? _rbp.getHGraph().getNodeHeight(i)
                              : _rbp.getHGraph().getNodeWidth(i);

      minx = std::min(minx, _rbp.getPlacement()[i].x);
      miny = std::min(miny, _rbp.getPlacement()[i].y);
      maxx = std::max(maxx, _rbp.getPlacement()[i].x + cellWidth);
      maxy = std::max(maxy, _rbp.getPlacement()[i].y + cellHeight);
    }
    unitsperx = (maxx - minx) / (double)_rbp.getParameters().xpixels;
    unitspery = (maxy - miny) / (double)_rbp.getParameters().ypixels;
    if (_rbp.getParameters().compressedPlotting) {
      horizPixelBuf = new vector<vector<bool> >(
          _rbp.getParameters().ypixels,
          vector<bool>(_rbp.getParameters().xpixels, false));
      vertPixelBuf = new vector<vector<bool> >(
          _rbp.getParameters().ypixels,
          vector<bool>(_rbp.getParameters().xpixels, false));
    }
    os << "# minx: " << minx << " maxx: " << maxx
       << " units per x pixel: " << unitsperx << endl
       << "# miny: " << miny << " maxy: " << maxy
       << " units per y pixel: " << unitspery << endl
       << "# optimized for resolution " << _rbp.getParameters().xpixels << "x"
       << _rbp.getParameters().ypixels << endl
       << endl;
  }
  unsigned compressedCells = 0;
  double maxCRHeight = _rbp.getMaxHeightCoreRow();
  for (j = 0; j < nodes.size(); j++) {
    i = nodes[j];
    const Orient& nodeOrient = _rbp.getPlacement().getOrient(i);
    unsigned angle = nodeOrient.getAngle();
    bool notRotated = angle == 0 || angle == 180;
    cellWidth = notRotated ? _rbp.getHGraph().getNodeWidth(i)
                           : _rbp.getHGraph().getNodeHeight(i);
    cellHeight = notRotated ? _rbp.getHGraph().getNodeHeight(i)
                            : _rbp.getHGraph().getNodeWidth(i);
    if (cellHeight <= maxCRHeight && _rbp.isFixed(i)) {
      haveAny = true;
      double offset = 0.0;
      if (cellWidth > cellHeight) {
        offset = 0.2 * cellHeight;
      } else {
        offset = 0.2 * cellWidth;
      }

      if (!_rbp.getParameters().adaptivePlotting ||
          !_rbp.getParameters().compressedPlotting)
        os << "# Cell Name " << _rbp.getHGraph().getNodeNameByIndex(i) << endl;

      if (fabs(_rbp.getPlacement()[i].x) < 1e+306 &&
          fabs(_rbp.getPlacement()[i].y) < 1e+306) {
        if (_rbp.getParameters().adaptivePlotting && cellWidth <= unitsperx &&
            cellHeight <= unitspery) {
          if (_rbp.getParameters().compressedPlotting) {
            unsigned xcoord = static_cast<unsigned>(
                ceil((_rbp.getPlacement()[i].x + 0.5 * cellWidth - minx) /
                     unitsperx));
            if (xcoord > 0) --xcoord;

            unsigned ycoord = static_cast<unsigned>(
                ceil((_rbp.getPlacement()[i].y + 0.5 * cellHeight - miny) /
                     unitspery));
            if (ycoord > 0) --ycoord;

            (*horizPixelBuf)[ycoord][xcoord] = true;
            ++compressedCells;
          } else {
            os << setw(10) << _rbp.getPlacement()[i].x << "  " << setw(10)
               << _rbp.getPlacement()[i].y + 0.5 * cellHeight << endl;
            os << setw(10) << _rbp.getPlacement()[i].x + cellWidth << "  "
               << setw(10) << _rbp.getPlacement()[i].y + 0.5 * cellHeight
               << endl
               << endl;
          }
        } else if (_rbp.getParameters().adaptivePlotting &&
                   cellWidth <= unitsperx) {
          if (_rbp.getParameters().compressedPlotting) {
            unsigned xcoord = static_cast<unsigned>(
                ceil((_rbp.getPlacement()[i].x + 0.5 * cellWidth - minx) /
                     unitsperx));
            if (xcoord > 0) --xcoord;

            double dbl_ymincoord =
                (_rbp.getPlacement()[i].y - miny) / unitspery;
            double dbl_ymincoord_floor = floor(dbl_ymincoord);
            unsigned ymincoord = static_cast<unsigned>(dbl_ymincoord_floor);
            if ((dbl_ymincoord - dbl_ymincoord_floor) > 0.9) ++ymincoord;

            double dbl_ymaxcoord =
                (_rbp.getPlacement()[i].y + cellHeight - miny) / unitspery;
            double dbl_ymaxcoord_ceil = ceil(dbl_ymaxcoord);
            unsigned ymaxcoord = static_cast<unsigned>(dbl_ymaxcoord_ceil) - 1;
            if ((dbl_ymaxcoord_ceil - dbl_ymaxcoord) > 0.9) --ymaxcoord;

            for (unsigned k = ymincoord; k <= ymaxcoord; ++k)
              (*vertPixelBuf)[k][xcoord] = true;
            ++compressedCells;
          } else {
            os << setw(10) << _rbp.getPlacement()[i].x + 0.5 * cellWidth << "  "
               << setw(10) << _rbp.getPlacement()[i].y << endl;
            os << setw(10) << _rbp.getPlacement()[i].x + 0.5 * cellWidth << "  "
               << setw(10) << _rbp.getPlacement()[i].y + cellHeight << endl
               << endl;
          }
        } else if (_rbp.getParameters().adaptivePlotting &&
                   cellHeight <= unitspery) {
          if (_rbp.getParameters().compressedPlotting) {
            double dbl_xmincoord =
                (_rbp.getPlacement()[i].x - minx) / unitsperx;
            double dbl_xmincoord_floor = floor(dbl_xmincoord);
            unsigned xmincoord = static_cast<unsigned>(dbl_xmincoord_floor);
            if ((dbl_xmincoord - dbl_xmincoord_floor) > 0.9) ++xmincoord;

            double dbl_xmaxcoord =
                (_rbp.getPlacement()[i].x + cellWidth - minx) / unitsperx;
            double dbl_xmaxcoord_ceil = ceil(dbl_xmaxcoord);
            unsigned xmaxcoord = static_cast<unsigned>(dbl_xmaxcoord_ceil) - 1;
            if ((dbl_xmaxcoord_ceil - dbl_xmaxcoord) > 0.9) --xmaxcoord;

            unsigned ycoord = static_cast<unsigned>(
                ceil((_rbp.getPlacement()[i].y + 0.5 * cellHeight - miny) /
                     unitspery));
            if (ycoord > 0) --ycoord;

            for (unsigned k = xmincoord; k <= xmaxcoord; ++k)
              (*horizPixelBuf)[ycoord][k] = true;
            ++compressedCells;
          } else {
            os << setw(10) << _rbp.getPlacement()[i].x << "  " << setw(10)
               << _rbp.getPlacement()[i].y + 0.5 * cellHeight << endl;
            os << setw(10) << _rbp.getPlacement()[i].x + cellWidth << "  "
               << setw(10) << _rbp.getPlacement()[i].y + 0.5 * cellHeight
               << endl
               << endl;
          }
        } else if (_rbp.getParameters().adaptivePlotting &&
                   (cellWidth <= 4. * unitsperx ||
                    cellHeight <= 4. * unitspery)) {
          if (_rbp.getParameters().compressedPlotting) {
            double dbl_xmincoord =
                (_rbp.getPlacement()[i].x - minx) / unitsperx;
            double dbl_xmincoord_floor = floor(dbl_xmincoord);
            unsigned xmincoord = static_cast<unsigned>(dbl_xmincoord_floor);
            if ((dbl_xmincoord - dbl_xmincoord_floor) > 0.9) ++xmincoord;

            double dbl_xmaxcoord =
                (_rbp.getPlacement()[i].x + cellWidth - minx) / unitsperx;
            double dbl_xmaxcoord_ceil = ceil(dbl_xmaxcoord);
            unsigned xmaxcoord = static_cast<unsigned>(dbl_xmaxcoord_ceil) - 1;
            if ((dbl_xmaxcoord_ceil - dbl_xmaxcoord) > 0.9) --xmaxcoord;

            double dbl_ymincoord =
                (_rbp.getPlacement()[i].y - miny) / unitspery;
            double dbl_ymincoord_floor = floor(dbl_ymincoord);
            unsigned ymincoord = static_cast<unsigned>(dbl_ymincoord_floor);
            if ((dbl_ymincoord - dbl_ymincoord_floor) > 0.9) ++ymincoord;

            double dbl_ymaxcoord =
                (_rbp.getPlacement()[i].y + cellHeight - miny) / unitspery;
            double dbl_ymaxcoord_ceil = ceil(dbl_ymaxcoord);
            unsigned ymaxcoord = static_cast<unsigned>(dbl_ymaxcoord_ceil) - 1;
            if ((dbl_ymaxcoord_ceil - dbl_ymaxcoord) > 0.9) --ymaxcoord;

            for (unsigned k = xmincoord; k <= xmaxcoord; ++k) {
              (*horizPixelBuf)[ymincoord][k] = true;
              (*horizPixelBuf)[ymaxcoord][k] = true;
            }
            for (unsigned k = ymincoord; k <= ymaxcoord; ++k) {
              (*vertPixelBuf)[k][xmincoord] = true;
              (*vertPixelBuf)[k][xmaxcoord] = true;
            }
            ++compressedCells;
          } else {
            os << setw(10) << _rbp.getPlacement()[i].x << "  " << setw(10)
               << _rbp.getPlacement()[i].y << endl;
            os << setw(10) << _rbp.getPlacement()[i].x + cellWidth << "  "
               << setw(10) << _rbp.getPlacement()[i].y << endl;
            os << setw(10) << _rbp.getPlacement()[i].x + cellWidth << "  "
               << setw(10) << _rbp.getPlacement()[i].y + cellHeight << endl;
            os << setw(10) << _rbp.getPlacement()[i].x << "  " << setw(10)
               << _rbp.getPlacement()[i].y + cellHeight << endl;
            os << setw(10) << _rbp.getPlacement()[i].x << "  " << setw(10)
               << _rbp.getPlacement()[i].y << endl
               << endl;
          }
        } else {
          if (nodeOrient == Orient("S")) {
            os << setw(10) << _rbp.getPlacement()[i].x << "  " << setw(10)
               << _rbp.getPlacement()[i].y + offset << endl;
            os << setw(10) << _rbp.getPlacement()[i].x + 3. * offset << "  "
               << setw(10) << _rbp.getPlacement()[i].y << endl;
          } else if (nodeOrient == Orient("FE")) {
            os << setw(10) << _rbp.getPlacement()[i].x << "  " << setw(10)
               << _rbp.getPlacement()[i].y + 3. * offset << endl;
            os << setw(10) << _rbp.getPlacement()[i].x + offset << "  "
               << setw(10) << _rbp.getPlacement()[i].y << endl;
          } else {
            os << setw(10) << _rbp.getPlacement()[i].x << "  " << setw(10)
               << _rbp.getPlacement()[i].y << endl;
          }

          if (nodeOrient == Orient("E")) {
            os << setw(10) << _rbp.getPlacement()[i].x + cellWidth - offset
               << "  " << setw(10) << _rbp.getPlacement()[i].y << endl;
            os << setw(10) << _rbp.getPlacement()[i].x + cellWidth << "  "
               << setw(10) << _rbp.getPlacement()[i].y + 3. * offset << endl;
          } else if (nodeOrient == Orient("FS")) {
            os << setw(10) << _rbp.getPlacement()[i].x + cellWidth - 3. * offset
               << "  " << setw(10) << _rbp.getPlacement()[i].y << endl;
            os << setw(10) << _rbp.getPlacement()[i].x + cellWidth << "  "
               << setw(10) << _rbp.getPlacement()[i].y + offset << endl;
          } else {
            os << setw(10) << _rbp.getPlacement()[i].x + cellWidth << "  "
               << setw(10) << _rbp.getPlacement()[i].y << endl;
          }

          if (nodeOrient == Orient("N")) {
            os << setw(10) << _rbp.getPlacement()[i].x + cellWidth << "  "
               << setw(10) << _rbp.getPlacement()[i].y + cellHeight - offset
               << endl;
            os << setw(10) << _rbp.getPlacement()[i].x + cellWidth - 3. * offset
               << "  " << setw(10) << _rbp.getPlacement()[i].y + cellHeight
               << endl;
          } else if (nodeOrient == Orient("FW")) {
            os << setw(10) << _rbp.getPlacement()[i].x + cellWidth << "  "
               << setw(10)
               << _rbp.getPlacement()[i].y + cellHeight - 3. * offset << endl;
            os << setw(10) << _rbp.getPlacement()[i].x + cellWidth - offset
               << "  " << setw(10) << _rbp.getPlacement()[i].y + cellHeight
               << endl;
          } else {
            os << setw(10) << _rbp.getPlacement()[i].x + cellWidth << "  "
               << setw(10) << _rbp.getPlacement()[i].y + cellHeight << endl;
          }

          if (nodeOrient == Orient("W")) {
            os << setw(10) << _rbp.getPlacement()[i].x + offset << "  "
               << setw(10) << _rbp.getPlacement()[i].y + cellHeight << endl;
            os << setw(10) << _rbp.getPlacement()[i].x << "  " << setw(10)
               << _rbp.getPlacement()[i].y + cellHeight - 3. * offset << endl;
          } else if (nodeOrient == Orient("FN")) {
            os << setw(10) << _rbp.getPlacement()[i].x + 3. * offset << "  "
               << setw(10) << _rbp.getPlacement()[i].y + cellHeight << endl;
            os << setw(10) << _rbp.getPlacement()[i].x << "  " << setw(10)
               << _rbp.getPlacement()[i].y + cellHeight - offset << endl;
          } else {
            os << setw(10) << _rbp.getPlacement()[i].x << "  " << setw(10)
               << _rbp.getPlacement()[i].y + cellHeight << endl;
          }

          if (nodeOrient == Orient("S")) {
            os << setw(10) << _rbp.getPlacement()[i].x << "  " << setw(10)
               << _rbp.getPlacement()[i].y + offset << endl
               << endl;
          } else if (nodeOrient == Orient("FE")) {
            os << setw(10) << _rbp.getPlacement()[i].x << "  " << setw(10)
               << _rbp.getPlacement()[i].y + 3. * offset << endl
               << endl;
          } else {
            os << setw(10) << _rbp.getPlacement()[i].x << "  " << setw(10)
               << _rbp.getPlacement()[i].y << endl
               << endl;
          }
        }
      }
    }
  }
  if (_rbp.getParameters().adaptivePlotting &&
      _rbp.getParameters().compressedPlotting) {
    os << "# Horizontal lines for " << compressedCells << " compressed cells"
       << endl
       << endl;
    for (i = 0; i < _rbp.getParameters().ypixels; ++i) {
      unsigned xmin = UINT_MAX;
      for (j = 0; j < _rbp.getParameters().xpixels; ++j) {
        if ((*horizPixelBuf)[i][j]) {
          (*vertPixelBuf)[i][j] = false;  // avoid drawing things twice
          if (xmin == UINT_MAX) xmin = j;
        } else {
          if (xmin != UINT_MAX) {
            os << setw(10) << minx + (double)xmin * unitsperx << "  "
               << setw(10) << miny + ((double)i + 0.5) * unitspery << endl;
            os << setw(10) << minx + (double)j * unitsperx << "  " << setw(10)
               << miny + ((double)i + 0.5) * unitspery << endl
               << endl;
            xmin = UINT_MAX;
          }
        }
      }
      if (xmin != UINT_MAX) {
        os << setw(10) << minx + (double)xmin * unitsperx << "  " << setw(10)
           << miny + ((double)i + 0.5) * unitspery << endl;
        os << setw(10)
           << minx + (double)_rbp.getParameters().xpixels * unitsperx << "  "
           << setw(10) << miny + ((double)i + 0.5) * unitspery << endl
           << endl;
      }
    }
    delete horizPixelBuf;
    os << "# Vertical lines for " << compressedCells << " compressed cells"
       << endl
       << endl;
    for (j = 0; j < _rbp.getParameters().xpixels; ++j) {
      unsigned ymin = UINT_MAX;
      for (i = 0; i < _rbp.getParameters().ypixels; ++i) {
        if ((*vertPixelBuf)[i][j]) {
          if (ymin == UINT_MAX) ymin = i;
        } else {
          if (ymin != UINT_MAX) {
            os << setw(10) << minx + ((double)j + 0.5) * unitsperx << "  "
               << setw(10) << miny + (double)ymin * unitspery << endl;
            os << setw(10) << minx + ((double)j + 0.5) * unitsperx << "  "
               << setw(10) << miny + (double)i * unitspery << endl
               << endl;
            ymin = UINT_MAX;
          }
        }
      }
      if (ymin != UINT_MAX) {
        os << setw(10) << minx + ((double)j + 0.5) * unitsperx << "  "
           << setw(10) << miny + (double)ymin * unitspery << endl;
        os << setw(10) << minx + ((double)j + 0.5) * unitsperx << "  "
           << setw(10)
           << miny + (double)_rbp.getParameters().ypixels * unitspery << endl
           << endl;
      }
    }
    delete vertPixelBuf;
  }

  if (!haveAny) os << "0" << endl << endl;

  os << "EOF" << endl << endl;
}

void RBPlacePlotter::writeGnuplotSmallMovableNodes(std::ostream& os,
                                                   const nodeSet& nodes) const {
  bool haveAny = false;
  vector<vector<bool> >*horizPixelBuf = 0, *vertPixelBuf = 0;
  double cellWidth = 0., cellHeight = 0.;
  unsigned i = 0, j;
  double unitsperx = 0.0, unitspery = 0.0;
  double minx = DBL_MAX, miny = DBL_MAX, maxx = -DBL_MAX, maxy = -DBL_MAX;
  if (_rbp.getParameters().adaptivePlotting) {
    for (j = 0; j < nodes.size(); j++) {
      i = nodes[j];
      if (!_modes.plotFillerCells &&
          _rbp.getHGraph().getNodeByIdx(i).getDegree() == 0)
        continue;
      if (fabs(_rbp.getPlacement()[i].x) >= 1e+306 ||
          fabs(_rbp.getPlacement()[i].y) >= 1e+306)
        continue;
      const Orient& nodeOrient = _rbp.getPlacement().getOrient(i);
      unsigned angle = nodeOrient.getAngle();
      bool notRotated = angle == 0 || angle == 180;
      cellWidth = notRotated ? _rbp.getHGraph().getNodeWidth(i)
                             : _rbp.getHGraph().getNodeHeight(i);
      cellHeight = notRotated ? _rbp.getHGraph().getNodeHeight(i)
                              : _rbp.getHGraph().getNodeWidth(i);
      minx = std::min(minx, _rbp.getPlacement()[i].x);
      miny = std::min(miny, _rbp.getPlacement()[i].y);
      maxx = std::max(maxx, _rbp.getPlacement()[i].x + cellWidth);
      maxy = std::max(maxy, _rbp.getPlacement()[i].y + cellHeight);
    }
    unitsperx = (maxx - minx) / (double)_rbp.getParameters().xpixels;
    unitspery = (maxy - miny) / (double)_rbp.getParameters().ypixels;
    if (_rbp.getParameters().compressedPlotting) {
      horizPixelBuf = new vector<vector<bool> >(
          _rbp.getParameters().ypixels,
          vector<bool>(_rbp.getParameters().xpixels, false));
      vertPixelBuf = new vector<vector<bool> >(
          _rbp.getParameters().ypixels,
          vector<bool>(_rbp.getParameters().xpixels, false));
    }
    os << "# minx: " << minx << " maxx: " << maxx
       << " units per x pixel: " << unitsperx << endl
       << "# miny: " << miny << " maxy: " << maxy
       << " units per y pixel: " << unitspery << endl
       << "# optimized for resolution " << _rbp.getParameters().xpixels << "x"
       << _rbp.getParameters().ypixels << endl
       << endl;
  }
  unsigned compressedCells = 0;
  double maxCRHeight = _rbp.getMaxHeightCoreRow();
  for (j = 0; j < nodes.size(); j++) {
    i = nodes[j];
    if (!_modes.plotFillerCells &&
        _rbp.getHGraph().getNodeByIdx(i).getDegree() == 0)
      continue;
    const Orient& nodeOrient = _rbp.getPlacement().getOrient(i);
    unsigned angle = nodeOrient.getAngle();
    bool notRotated = angle == 0 || angle == 180;
    cellWidth = notRotated ? _rbp.getHGraph().getNodeWidth(i)
                           : _rbp.getHGraph().getNodeHeight(i);
    cellHeight = notRotated ? _rbp.getHGraph().getNodeHeight(i)
                            : _rbp.getHGraph().getNodeWidth(i);
    if (cellHeight <= maxCRHeight && !_rbp.isFixed(i)) {
      haveAny = true;
      double offset = 0.0;
      if (cellWidth > cellHeight) {
        offset = 0.2 * cellHeight;
      } else {
        offset = 0.2 * cellWidth;
      }

      if (!_rbp.getParameters().adaptivePlotting ||
          !_rbp.getParameters().compressedPlotting)
        os << "# Cell Name " << _rbp.getHGraph().getNodeNameByIndex(i) << endl;

      if (fabs(_rbp.getPlacement()[i].x) < 1e+306 &&
          fabs(_rbp.getPlacement()[i].y) < 1e+306) {
        if (_rbp.getParameters().adaptivePlotting && cellWidth <= unitsperx &&
            cellHeight <= unitspery) {
          if (_rbp.getParameters().compressedPlotting) {
            unsigned xcoord = static_cast<unsigned>(
                ceil((_rbp.getPlacement()[i].x + 0.5 * cellWidth - minx) /
                     unitsperx));
            if (xcoord > 0) --xcoord;

            unsigned ycoord = static_cast<unsigned>(
                ceil((_rbp.getPlacement()[i].y + 0.5 * cellHeight - miny) /
                     unitspery));
            if (ycoord > 0) --ycoord;

            (*horizPixelBuf)[ycoord][xcoord] = true;
            ++compressedCells;
          } else {
            os << setw(10) << _rbp.getPlacement()[i].x << "  " << setw(10)
               << _rbp.getPlacement()[i].y + 0.5 * cellHeight << endl;
            os << setw(10) << _rbp.getPlacement()[i].x + cellWidth << "  "
               << setw(10) << _rbp.getPlacement()[i].y + 0.5 * cellHeight
               << endl
               << endl;
          }
        } else if (_rbp.getParameters().adaptivePlotting &&
                   cellWidth <= unitsperx) {
          if (_rbp.getParameters().compressedPlotting) {
            unsigned xcoord = static_cast<unsigned>(
                ceil((_rbp.getPlacement()[i].x + 0.5 * cellWidth - minx) /
                     unitsperx));
            if (xcoord > 0) --xcoord;

            double dbl_ymincoord =
                (_rbp.getPlacement()[i].y - miny) / unitspery;
            double dbl_ymincoord_floor = floor(dbl_ymincoord);
            unsigned ymincoord = static_cast<unsigned>(dbl_ymincoord_floor);
            if ((dbl_ymincoord - dbl_ymincoord_floor) > 0.9) ++ymincoord;

            double dbl_ymaxcoord =
                (_rbp.getPlacement()[i].y + cellHeight - miny) / unitspery;
            double dbl_ymaxcoord_ceil = ceil(dbl_ymaxcoord);
            unsigned ymaxcoord = static_cast<unsigned>(dbl_ymaxcoord_ceil) - 1;
            if ((dbl_ymaxcoord_ceil - dbl_ymaxcoord) > 0.9) --ymaxcoord;

            for (unsigned k = ymincoord; k <= ymaxcoord; ++k)
              (*vertPixelBuf)[k][xcoord] = true;
            ++compressedCells;
          } else {
            os << setw(10) << _rbp.getPlacement()[i].x + 0.5 * cellWidth << "  "
               << setw(10) << _rbp.getPlacement()[i].y << endl;
            os << setw(10) << _rbp.getPlacement()[i].x + 0.5 * cellWidth << "  "
               << setw(10) << _rbp.getPlacement()[i].y + cellHeight << endl
               << endl;
          }
        } else if (_rbp.getParameters().adaptivePlotting &&
                   cellHeight <= unitspery) {
          if (_rbp.getParameters().compressedPlotting) {
            double dbl_xmincoord =
                (_rbp.getPlacement()[i].x - minx) / unitsperx;
            double dbl_xmincoord_floor = floor(dbl_xmincoord);
            unsigned xmincoord = static_cast<unsigned>(dbl_xmincoord_floor);
            if ((dbl_xmincoord - dbl_xmincoord_floor) > 0.9) ++xmincoord;

            double dbl_xmaxcoord =
                (_rbp.getPlacement()[i].x + cellWidth - minx) / unitsperx;
            double dbl_xmaxcoord_ceil = ceil(dbl_xmaxcoord);
            unsigned xmaxcoord = static_cast<unsigned>(dbl_xmaxcoord_ceil) - 1;
            if ((dbl_xmaxcoord_ceil - dbl_xmaxcoord) > 0.9) --xmaxcoord;

            unsigned ycoord = static_cast<unsigned>(
                ceil((_rbp.getPlacement()[i].y + 0.5 * cellHeight - miny) /
                     unitspery));
            if (ycoord > 0) --ycoord;

            for (unsigned k = xmincoord; k <= xmaxcoord; ++k)
              (*horizPixelBuf)[ycoord][k] = true;
            ++compressedCells;
          } else {
            os << setw(10) << _rbp.getPlacement()[i].x << "  " << setw(10)
               << _rbp.getPlacement()[i].y + 0.5 * cellHeight << endl;
            os << setw(10) << _rbp.getPlacement()[i].x + cellWidth << "  "
               << setw(10) << _rbp.getPlacement()[i].y + 0.5 * cellHeight
               << endl
               << endl;
          }
        } else if (_rbp.getParameters().adaptivePlotting &&
                   (cellWidth <= 4. * unitsperx ||
                    cellHeight <= 4. * unitspery)) {
          if (_rbp.getParameters().compressedPlotting) {
            double dbl_xmincoord =
                (_rbp.getPlacement()[i].x - minx) / unitsperx;
            double dbl_xmincoord_floor = floor(dbl_xmincoord);
            unsigned xmincoord = static_cast<unsigned>(dbl_xmincoord_floor);
            if ((dbl_xmincoord - dbl_xmincoord_floor) > 0.9) ++xmincoord;

            double dbl_xmaxcoord =
                (_rbp.getPlacement()[i].x + cellWidth - minx) / unitsperx;
            double dbl_xmaxcoord_ceil = ceil(dbl_xmaxcoord);
            unsigned xmaxcoord = static_cast<unsigned>(dbl_xmaxcoord_ceil) - 1;
            if ((dbl_xmaxcoord_ceil - dbl_xmaxcoord) > 0.9) --xmaxcoord;

            double dbl_ymincoord =
                (_rbp.getPlacement()[i].y - miny) / unitspery;
            double dbl_ymincoord_floor = floor(dbl_ymincoord);
            unsigned ymincoord = static_cast<unsigned>(dbl_ymincoord_floor);
            if ((dbl_ymincoord - dbl_ymincoord_floor) > 0.9) ++ymincoord;

            double dbl_ymaxcoord =
                (_rbp.getPlacement()[i].y + cellHeight - miny) / unitspery;
            double dbl_ymaxcoord_ceil = ceil(dbl_ymaxcoord);
            unsigned ymaxcoord = static_cast<unsigned>(dbl_ymaxcoord_ceil) - 1;
            if ((dbl_ymaxcoord_ceil - dbl_ymaxcoord) > 0.9) --ymaxcoord;

            for (unsigned k = xmincoord; k <= xmaxcoord; ++k) {
              (*horizPixelBuf)[ymincoord][k] = true;
              (*horizPixelBuf)[ymaxcoord][k] = true;
            }
            for (unsigned k = ymincoord; k <= ymaxcoord; ++k) {
              (*vertPixelBuf)[k][xmincoord] = true;
              (*vertPixelBuf)[k][xmaxcoord] = true;
            }
            ++compressedCells;
          } else {
            os << setw(10) << _rbp.getPlacement()[i].x << "  " << setw(10)
               << _rbp.getPlacement()[i].y << endl;
            os << setw(10) << _rbp.getPlacement()[i].x + cellWidth << "  "
               << setw(10) << _rbp.getPlacement()[i].y << endl;
            os << setw(10) << _rbp.getPlacement()[i].x + cellWidth << "  "
               << setw(10) << _rbp.getPlacement()[i].y + cellHeight << endl;
            os << setw(10) << _rbp.getPlacement()[i].x << "  " << setw(10)
               << _rbp.getPlacement()[i].y + cellHeight << endl;
            os << setw(10) << _rbp.getPlacement()[i].x << "  " << setw(10)
               << _rbp.getPlacement()[i].y << endl
               << endl;
          }
        } else {
          if (nodeOrient == Orient("S")) {
            os << setw(10) << _rbp.getPlacement()[i].x << "  " << setw(10)
               << _rbp.getPlacement()[i].y + offset << endl;
            os << setw(10) << _rbp.getPlacement()[i].x + 3. * offset << "  "
               << setw(10) << _rbp.getPlacement()[i].y << endl;
          } else if (nodeOrient == Orient("FE")) {
            os << setw(10) << _rbp.getPlacement()[i].x << "  " << setw(10)
               << _rbp.getPlacement()[i].y + 3. * offset << endl;
            os << setw(10) << _rbp.getPlacement()[i].x + offset << "  "
               << setw(10) << _rbp.getPlacement()[i].y << endl;
          } else {
            os << setw(10) << _rbp.getPlacement()[i].x << "  " << setw(10)
               << _rbp.getPlacement()[i].y << endl;
          }

          if (nodeOrient == Orient("E")) {
            os << setw(10) << _rbp.getPlacement()[i].x + cellWidth - offset
               << "  " << setw(10) << _rbp.getPlacement()[i].y << endl;
            os << setw(10) << _rbp.getPlacement()[i].x + cellWidth << "  "
               << setw(10) << _rbp.getPlacement()[i].y + 3. * offset << endl;
          } else if (nodeOrient == Orient("FS")) {
            os << setw(10) << _rbp.getPlacement()[i].x + cellWidth - 3. * offset
               << "  " << setw(10) << _rbp.getPlacement()[i].y << endl;
            os << setw(10) << _rbp.getPlacement()[i].x + cellWidth << "  "
               << setw(10) << _rbp.getPlacement()[i].y + offset << endl;
          } else {
            os << setw(10) << _rbp.getPlacement()[i].x + cellWidth << "  "
               << setw(10) << _rbp.getPlacement()[i].y << endl;
          }

          if (nodeOrient == Orient("N")) {
            os << setw(10) << _rbp.getPlacement()[i].x + cellWidth << "  "
               << setw(10) << _rbp.getPlacement()[i].y + cellHeight - offset
               << endl;
            os << setw(10) << _rbp.getPlacement()[i].x + cellWidth - 3. * offset
               << "  " << setw(10) << _rbp.getPlacement()[i].y + cellHeight
               << endl;
          } else if (nodeOrient == Orient("FW")) {
            os << setw(10) << _rbp.getPlacement()[i].x + cellWidth << "  "
               << setw(10)
               << _rbp.getPlacement()[i].y + cellHeight - 3. * offset << endl;
            os << setw(10) << _rbp.getPlacement()[i].x + cellWidth - offset
               << "  " << setw(10) << _rbp.getPlacement()[i].y + cellHeight
               << endl;
          } else {
            os << setw(10) << _rbp.getPlacement()[i].x + cellWidth << "  "
               << setw(10) << _rbp.getPlacement()[i].y + cellHeight << endl;
          }

          if (nodeOrient == Orient("W")) {
            os << setw(10) << _rbp.getPlacement()[i].x + offset << "  "
               << setw(10) << _rbp.getPlacement()[i].y + cellHeight << endl;
            os << setw(10) << _rbp.getPlacement()[i].x << "  " << setw(10)
               << _rbp.getPlacement()[i].y + cellHeight - 3. * offset << endl;
          } else if (nodeOrient == Orient("FN")) {
            os << setw(10) << _rbp.getPlacement()[i].x + 3. * offset << "  "
               << setw(10) << _rbp.getPlacement()[i].y + cellHeight << endl;
            os << setw(10) << _rbp.getPlacement()[i].x << "  " << setw(10)
               << _rbp.getPlacement()[i].y + cellHeight - offset << endl;
          } else {
            os << setw(10) << _rbp.getPlacement()[i].x << "  " << setw(10)
               << _rbp.getPlacement()[i].y + cellHeight << endl;
          }

          if (nodeOrient == Orient("S")) {
            os << setw(10) << _rbp.getPlacement()[i].x << "  " << setw(10)
               << _rbp.getPlacement()[i].y + offset << endl
               << endl;
          } else if (nodeOrient == Orient("FE")) {
            os << setw(10) << _rbp.getPlacement()[i].x << "  " << setw(10)
               << _rbp.getPlacement()[i].y + 3. * offset << endl
               << endl;
          } else {
            os << setw(10) << _rbp.getPlacement()[i].x << "  " << setw(10)
               << _rbp.getPlacement()[i].y << endl
               << endl;
          }
        }
      }
    }
  }
  if (_rbp.getParameters().adaptivePlotting &&
      _rbp.getParameters().compressedPlotting) {
    os << "# Horizontal lines for " << compressedCells << " compressed cells"
       << endl
       << endl;
    for (i = 0; i < _rbp.getParameters().ypixels; ++i) {
      unsigned xmin = UINT_MAX;
      for (j = 0; j < _rbp.getParameters().xpixels; ++j) {
        if ((*horizPixelBuf)[i][j]) {
          (*vertPixelBuf)[i][j] = false;  // avoid drawing things twice
          if (xmin == UINT_MAX) xmin = j;
        } else {
          if (xmin != UINT_MAX) {
            os << setw(10) << minx + (double)xmin * unitsperx << "  "
               << setw(10) << miny + ((double)i + 0.5) * unitspery << endl;
            os << setw(10) << minx + (double)j * unitsperx << "  " << setw(10)
               << miny + ((double)i + 0.5) * unitspery << endl
               << endl;
            xmin = UINT_MAX;
          }
        }
      }
      if (xmin != UINT_MAX) {
        os << setw(10) << minx + (double)xmin * unitsperx << "  " << setw(10)
           << miny + ((double)i + 0.5) * unitspery << endl;
        os << setw(10)
           << minx + (double)_rbp.getParameters().xpixels * unitsperx << "  "
           << setw(10) << miny + ((double)i + 0.5) * unitspery << endl
           << endl;
      }
    }
    delete horizPixelBuf;
    os << "# Vertical lines for " << compressedCells << " compressed cells"
       << endl
       << endl;
    for (j = 0; j < _rbp.getParameters().xpixels; ++j) {
      unsigned ymin = UINT_MAX;
      for (i = 0; i < _rbp.getParameters().ypixels; ++i) {
        if ((*vertPixelBuf)[i][j]) {
          if (ymin == UINT_MAX) ymin = i;
        } else {
          if (ymin != UINT_MAX) {
            os << setw(10) << minx + ((double)j + 0.5) * unitsperx << "  "
               << setw(10) << miny + (double)ymin * unitspery << endl;
            os << setw(10) << minx + ((double)j + 0.5) * unitsperx << "  "
               << setw(10) << miny + (double)i * unitspery << endl
               << endl;
            ymin = UINT_MAX;
          }
        }
      }
      if (ymin != UINT_MAX) {
        os << setw(10) << minx + ((double)j + 0.5) * unitsperx << "  "
           << setw(10) << miny + (double)ymin * unitspery << endl;
        os << setw(10) << minx + ((double)j + 0.5) * unitsperx << "  "
           << setw(10)
           << miny + (double)_rbp.getParameters().ypixels * unitspery << endl
           << endl;
      }
    }
    delete vertPixelBuf;
  }

  if (!haveAny) os << "0" << endl << endl;

  os << "EOF" << endl << endl;
}

void RBPlacePlotter::writeGnuplotSiteMap(std::ostream& os) const {
  double rowYCoord = 0;
  double siteWidth = 0;
  double rowHeight = 0;
  double rowSpacing = 0;
  unsigned numSitesSubRow = 0;

  itRBPCoreRow CR;
  itRBPSubRow SR;
  for (CR = _rbp.coreRowsBegin(); CR != _rbp.coreRowsEnd(); ++CR) {
    rowYCoord = (CR)->getYCoord();
    siteWidth = (CR)->getSiteWidth();
    rowHeight = (CR)->getHeight();
    rowSpacing = (CR)->getSpacing();

    for (SR = (CR)->subRowsBegin(); SR != (CR)->subRowsEnd(); ++SR) {
      double xStart = SR->getXStart();
      numSitesSubRow = (SR)->getNumSites();
      unsigned i = 0;
      for (i = 0; i < numSitesSubRow; ++i) {
        double siteXstart = xStart + i * rowSpacing;
        double siteXend = siteXstart + siteWidth;

        os << setw(10) << siteXstart << "  " << (rowYCoord + 0.1 * rowHeight)
           << endl;
        os << setw(10) << siteXstart + 0.1 * siteWidth << "  " << (rowYCoord)
           << endl;
        os << setw(10) << siteXend - 0.1 * siteWidth << "  " << (rowYCoord)
           << endl;
        os << setw(10) << siteXend << "  " << (rowYCoord + 0.1 * rowHeight)
           << endl
           << endl;
      }
    }
  }
  os << "EOF" << endl;
}

void RBPlacePlotter::writeGnuplotRows(std::ostream& os) const {
  double rowYCoord = 0;
  double siteWidth = 0;
  double rowHeight = 0;
  unsigned numSitesSubRow = 0;

  itRBPCoreRow CR;
  itRBPSubRow SR;
  for (CR = _rbp.coreRowsBegin(); CR != _rbp.coreRowsEnd(); ++CR) {
    rowYCoord = (CR)->getYCoord();
    siteWidth = (CR)->getSiteWidth();
    rowHeight = (CR)->getHeight();
    for (SR = (CR)->subRowsBegin(); SR != (CR)->subRowsEnd(); ++SR) {
      double xStart = SR->getXStart();
      numSitesSubRow = (SR)->getNumSites();
      double xEnd = xStart + numSitesSubRow * siteWidth;

      os << setw(10) << xStart << "  " << (rowYCoord + 0.1 * rowHeight) << endl;
      os << setw(10) << xStart + 0.1 * siteWidth << "  " << (rowYCoord) << endl;
      os << setw(10) << xEnd - 0.1 * siteWidth << "  " << (rowYCoord) << endl;
      os << setw(10) << xEnd << "  " << (rowYCoord + 0.1 * rowHeight) << endl
         << endl;
    }
  }
  os << "EOF" << endl;
}

size_t RBPlacePlotter::getFileCount(void) const {
  size_t rval = 0;
  if (_modes.plotNodes) rval += _nodeGroups.size() * 4;
  if (_modes.plotNets) rval += _nodes4Nets.size();
  if (_modes.plotSiteMap) rval += 1;
  if (_modes.plotRows) rval += 1;

  _colors.reserve(rval);

  return rval;
}

void RBPlacePlotter::addNodeSet(const nodeSet& ns) {
  _nodeGroups.push_back(ns);
  update();
}

void RBPlacePlotter::addNodeSetForLabels(const nodeSet& ns) {
  _nodes2Label.push_back(ns);
  update();
}

void RBPlacePlotter::addNodeSetForNets(const nodeSet& ns) {
  _nodes4Nets.push_back(ns);
  update();
}

void RBPlacePlotter::assignRandomColors(void) const {
  RandomRawUnsigned r;
  _colors.clear();
  for (size_t it = 0; it < _file_ct; ++it)
    _colors.push_back(Plotters::Colors(r % 7 + 1));
}

void RBPlacePlotter::assignStandardColors(void) const {
  _colors.clear();
  if (_modes.oneColorPerGroup) {
    _colors.push_back(Plotters::orange);  // Small Movable
    _colors.push_back(Plotters::orange);  // Small Fixed
    _colors.push_back(Plotters::orange);  // Large Movable
    _colors.push_back(Plotters::orange);  // Large Fixed
    Plotters::Colors color = Plotters::red;
    for (unsigned it = 1; it < _nodeGroups.size(); ++it) {
      _colors.push_back(color);  // Small Movable
      _colors.push_back(color);  // Small Fixed
      _colors.push_back(color);  // Large Movable
      _colors.push_back(color);  // Large Fixed
      color = static_cast<Plotters::Colors>(
          (static_cast<unsigned>(color) + 1u) % 9u);
    }
  } else {
    for (unsigned it = 0; it < _nodeGroups.size(); ++it) {
      _colors.push_back(Plotters::blue);    // Small Movable
      _colors.push_back(Plotters::orange);  // Small Fixed
      _colors.push_back(Plotters::green);   // Large Movable
      _colors.push_back(Plotters::orange);  // Large Fixed
    }
  }
  for (unsigned it = 0; it < _nodes4Nets.size(); ++it) {
    if (_modes.plotSteinerNets) {
      _colors.push_back(Plotters::red);
      _colors.push_back(Plotters::cyan);
      _colors.push_back(Plotters::sienna);
    } else {
      _colors.push_back(Plotters::blackdotted);
    }
  }
  if (_modes.plotSiteMap) _colors.push_back(Plotters::magenta);
  if (_modes.plotRows) _colors.push_back(Plotters::magenta);
}

}  // namespace Plotters

// END RBPLACEPLOTTER CODE

int set_of(unsigned idx, vector<int> tree) {
  int i = -90;
  while ((i = tree[idx]) != -1) idx = tree[idx];
  return idx;
}

void join(unsigned idx1, unsigned idx2, vector<int>& tree, bool coin) {
  int set1 = set_of(idx1, tree), set2 = set_of(idx2, tree);
  abkassert(tree[set1] == -1,
            "Invalid tree look up, Internal RBPlacePlotter Error");
  abkassert(tree[set2] == -1,
            "Invalid tree look up, Internal RBPlacePlotter Error");
  if (set1 != set2) {
    if (coin)
      tree[set1] = set2;
    else
      tree[set2] = set1;
  }
}

// The following function is only for debugging.
// It prints a "tree" which is stored in a vector
// It is intended to be used with my implementation of
// the algorithm to find connected components in a graph //TODO CITE SEDGEWICK
// ALGORITHMS PART 1-4, Ch. 1
void spam_tree(const vector<int>& tree) {
  cout << "Tree: ";
  for (size_t i = 0; i < tree.size(); ++i) cout << setw(10) << tree[i];
  cout << endl;
}
// BEGIN RBPlacePlotter consumer code
void Plotters::RBPlacePlotter::plotCellGroups(string baseFileName,
                                              const RBPlacement& rbp,
                                              unsigned numCellGroups) {
  RBPlacePlotter::Parameters plParams;
  plParams.plotNodes = true;
  plParams.oneColorPerGroup = true;
  RBPlacePlotter plot(plParams, baseFileName, rbp);

  const HGraphWDimensions& hGraph = rbp.getHGraph();
  vector<vector<unsigned> > cellGroups(numCellGroups + 1);
  for (unsigned i = 0; i < hGraph.getNumNodes(); ++i) {
    if (i < hGraph.getNumTerminals()) {
      cellGroups[0].push_back(i);
    } else {
      // determine the group number of the cell
      std::string cellName = hGraph.getNodeNameByIndex(i).c_str();
      std::string::size_type sepPos = cellName.rfind('_');
      std::string groupNumString = cellName.substr(sepPos + 1);
      // strtol (not atoi) so a non-numeric suffix such as "5x" is rejected
      // rather than silently truncated.
      errno = 0;
      char* gEnd = 0;
      long gid = std::strtol(groupNumString.c_str(), &gEnd, 10);
      abkfatal(gEnd != groupNumString.c_str() && *gEnd == '\0' && errno == 0,
               "non-numeric cell-group suffix in cell name");
      unsigned groupId = static_cast<unsigned>(gid);
      abkfatal(groupId <= numCellGroups, "bad groupid");
      abkfatal(groupId > 0, "bad groupid");
      cellGroups[groupId].push_back(i);
    }
  }

  for (unsigned i = 0; i <= numCellGroups; ++i) {
    plot.addNodeSet(cellGroups[i]);
  }

  ofstream gpFile;
  gpFile.open((string(baseFileName) + string(".gpl")).c_str());
  gpFile << plot;
  gpFile.close();
}

void Plotters::RBPlacePlotter::plotFreeshapeFloorplan(string baseFileName,
                                                      const RBPlacement& rbp,
                                                      int argc,
                                                      const char** argv) {
  RBPlacePlotter::Parameters gplParams(argc, argv);
  gplParams.plotNodes = true;
  gplParams.plotNets = false;
  // gplParams.plotNets=true;
  RBPlacePlotter plot(gplParams, baseFileName, rbp);
  RandomRawUnsigned r;

  const HGraphWDimensions& hGraph = rbp.getHGraph();
  vector<bool> v_visited(rbp.getHGraph().getNumNodes());
  // for all nets not visited n
  // seed nodeset with unvisited node from n
  // while(nodeset.size() == last_size)
  //   add everything unvisited connected via fake wire to every node unvisited

  vector<int> tree(hGraph.getNumNodes(), -1);
  std::set<unsigned> setForNets;
  for (itHGFEdgeGlobal it = hGraph.edgesBegin(); it != hGraph.edgesEnd();
       ++it) {
    string netName(hGraph.getNetNameByIndex((*it)->getIndex()));
    if (netName.find("temp") == string::npos) continue;
    const HGFEdge& this_net = *(*it);
    itHGFNodeLocal firstIt = this_net.nodesBegin();
    const HGFNode& firstNode = *(*firstIt);
    unsigned firstIdx = firstNode.getIndex();

    setForNets.insert(firstIdx);
    for (itHGFNodeLocal nodeIt = (*it)->nodesBegin();
         nodeIt != (*it)->nodesEnd(); ++nodeIt) {
      unsigned v = (*nodeIt)->getIndex();
      if (firstIdx != v) join(firstIdx, v, tree, r & 1u);
      setForNets.insert(v);
      // spam_tree(tree);
    }
  }

  //    Commented out to turn off nets plotting
  //    Plotters::nodeSet nsForNet;
  //    for(set<unsigned>::iterator it =
  //    setForNets.begin();it!=setForNets.end();++it)
  //    {
  //       nsForNet.push_back(*it);
  //    }
  //    plot.addNodeSetForNets(nsForNet);

  map<unsigned, Plotters::nodeSet> nodeSets;
  for (size_t node = 0; node != tree.size(); ++node) {
    unsigned node_set_num = set_of(node, tree);
    nodeSets[node_set_num].push_back(node);
  }

  // aggregate all single node sets into one large set
  Plotters::nodeSet singleNodes;
  std::set<unsigned> delete_us;
  for (map<unsigned, nodeSet>::const_iterator it = nodeSets.begin();
       it != nodeSets.end(); ++it)
    if (it->second.size() == 1) {
      singleNodes.push_back((it->second)[0]);
      delete_us.insert(it->first);
    }
  for (set<unsigned>::iterator it = delete_us.begin(); it != delete_us.end();
       ++it) {
    nodeSets.erase(*it);
  }
  if (singleNodes.size() > 0) plot.addNodeSet(singleNodes);
  // end aggregation

  for (map<unsigned, nodeSet>::const_iterator it = nodeSets.begin();
       it != nodeSets.end(); ++it) {
    if (it->second.size() == 0) abkfatal(0, "Cannot plot an empty node set\n");
    plot.addNodeSet(it->second);
  }

  ofstream gpFile;
  gpFile.open((string(baseFileName) + string(".gpl")).c_str());
  gpFile << plot;
  gpFile.close();
}
// END RBPlacePlotter consumer code

// BEGIN OLDER FUNCTIONS
// these are not necessarily depricated and ported to the above framework yet
struct ctrSortWSMap {
  unsigned horizIdx;
  unsigned vertIdx;
  double WS;
};

struct wsMap_less_mag {
  bool operator()(const ctrSortWSMap& elem1, const ctrSortWSMap& elem2) {
    return (elem1.WS < elem2.WS);
  }
};

void RBPlacement::plotWSHist(const char* baseFileName) {
  cout << "Saving " << baseFileName << ".gpl" << endl;

  string fname = string(baseFileName) + ".gpl";
  ofstream gpFile(fname.c_str());

  const HGraphWDimensions& hgWDims = getHGraph();
  unsigned numCells = hgWDims.getNumNodes() - hgWDims.getNumTerminals();
  BBox layoutBBox = getBBox(false);
  double xSize = layoutBBox.xMax - layoutBBox.xMin;
  double ySize = layoutBBox.yMax - layoutBBox.yMin;
  double horizGridSize =
      10 * sqrt(xSize * ySize / numCells);  // try to have 100 cells
  double vertGridSize = horizGridSize;      // in each grid box

  double rowSpacing =
      (coreRowsBegin() + 1)->getYCoord() - coreRowsBegin()->getYCoord();
  horizGridSize = ceil(horizGridSize / rowSpacing) * rowSpacing;
  vertGridSize = horizGridSize;

  if (horizGridSize > xSize) horizGridSize = xSize / 2;
  if (vertGridSize > ySize) vertGridSize = ySize / 2;

  unsigned maxHorizIdx = unsigned(ceil(xSize / horizGridSize));
  unsigned maxVertIdx = unsigned(ceil(ySize / vertGridSize));

  vector<vector<double> > wsMap(maxVertIdx);
  for (unsigned i = 0; i < maxVertIdx; ++i) {
    wsMap[i].resize(maxHorizIdx, 0);
  }
  vector<ctrSortWSMap> sortWSBins;
  sortWSBins.resize(maxVertIdx * maxHorizIdx);

  for (unsigned i = hgWDims.getNumTerminals(); i < hgWDims.getNumNodes(); ++i) {
    //    const HGFNode& node=hgWDims.getNodeByIdx(i);
    BBox nodeBBox;
    nodeBBox.xMin = _placement[i].x;
    nodeBBox.yMin = _placement[i].y;
    nodeBBox.xMax = nodeBBox.xMin + hgWDims.getNodeWidth(i);
    nodeBBox.yMax = nodeBBox.yMin + hgWDims.getNodeHeight(i);

    if (nodeBBox.xMin < layoutBBox.xMin) nodeBBox.xMin = layoutBBox.xMin;
    if (nodeBBox.xMax > layoutBBox.xMax) nodeBBox.xMax = layoutBBox.xMax;
    if (nodeBBox.yMin < layoutBBox.yMin) nodeBBox.yMin = layoutBBox.yMin;
    if (nodeBBox.yMax > layoutBBox.yMax) nodeBBox.yMax = layoutBBox.yMax;
    if (nodeBBox.xMax < layoutBBox.xMin) nodeBBox.xMax = layoutBBox.xMin;
    if (nodeBBox.xMin > layoutBBox.xMax) nodeBBox.xMin = layoutBBox.xMax;
    if (nodeBBox.yMax < layoutBBox.yMin) nodeBBox.yMax = layoutBBox.yMin;
    if (nodeBBox.yMin > layoutBBox.yMax) nodeBBox.yMin = layoutBBox.yMax;

    int minHorizIndex =
        int(floor((nodeBBox.xMin - layoutBBox.xMin) / horizGridSize));
    int maxHorizIndex =
        int(ceil((nodeBBox.xMax - layoutBBox.xMin) / horizGridSize));
    int minVertIndex =
        int(floor((nodeBBox.yMin - layoutBBox.yMin) / vertGridSize));
    int maxVertIndex =
        int(ceil((nodeBBox.yMax - layoutBBox.yMin) / vertGridSize));

    for (int k = minVertIndex; k < maxVertIndex; ++k) {
      for (int j = minHorizIndex; j < maxHorizIndex; ++j) {
        double regXMin = j * horizGridSize + layoutBBox.xMin;
        double regXMax = (j + 1) * horizGridSize + layoutBBox.xMin;
        double regYMin = k * vertGridSize + layoutBBox.yMin;
        double regYMax = (k + 1) * vertGridSize + layoutBBox.yMin;

        if (regXMax > layoutBBox.xMax) regXMax = layoutBBox.xMax;
        if (regYMax > layoutBBox.yMax) regYMax = layoutBBox.yMax;

        BBox region;
        region.add(regXMin, regYMin);
        region.add(regXMax, regYMax);

        BBox intersect = nodeBBox.intersect(region);
        double area = intersect.getHeight() * intersect.getWidth();
        wsMap[k][j] += area;
      }
    }
  }

  for (unsigned k = 0; k < maxVertIdx; ++k) {
    for (unsigned j = 0; j < maxHorizIdx; ++j) {
      double regXMin = j * horizGridSize + layoutBBox.xMin;
      double regXMax = (j + 1) * horizGridSize + layoutBBox.xMin;
      double regYMin = k * vertGridSize + layoutBBox.yMin;
      double regYMax = (k + 1) * vertGridSize + layoutBBox.yMin;

      if (regXMax > layoutBBox.xMax) regXMax = layoutBBox.xMax;
      if (regYMax > layoutBBox.yMax) regYMax = layoutBBox.yMax;

      BBox region;
      region.add(regXMin, regYMin);
      region.add(regXMax, regYMax);

      // double regionArea = region.getArea();
      double regionArea = getContainedSiteAreaInBBox(region);
      double regionWS = 1;
      if (regionArea != 0) regionWS = (regionArea - wsMap[k][j]) / regionArea;
      wsMap[k][j] = regionWS;
    }
  }

  for (unsigned i = 0; i < maxVertIdx; ++i) {
    for (unsigned j = 0; j < maxHorizIdx; ++j) {
      unsigned idx = i * maxHorizIdx + j;
      sortWSBins[idx].horizIdx = j;
      sortWSBins[idx].vertIdx = i;
      sortWSBins[idx].WS = wsMap[i][j];
    }
  }
  sort(sortWSBins.begin(), sortWSBins.end(), wsMap_less_mag());

  int numWSBins = maxVertIdx * maxHorizIdx;

  cout << "Whitespace Map Statistics " << endl;
  cout << "Core Layout BBox " << layoutBBox;
  cout << "Whitespace grid is " << maxVertIdx << " x " << maxHorizIdx << " = "
       << numWSBins << " bins" << endl;
  double maxWS = sortWSBins[numWSBins - 1].WS;
  double minWS = sortWSBins[0].WS;
  double avgWS = 0;
  double medianWS = sortWSBins[int(numWSBins / 2)].WS;

  for (int i = 0; i < numWSBins; ++i) avgWS += sortWSBins[i].WS;
  avgWS /= numWSBins;

  int numLevels = 20;
  double rangeLevel = 1.0 / double(numLevels);
  vector<int> binsPerLevel(numLevels, 0);

  double currLevel = 0 + rangeLevel;
  int currLevelIdx = 0;
  int currIdx = 0;
  while (currIdx < numWSBins) {
    while (currIdx < numWSBins) {
      if (sortWSBins[currIdx].WS <= currLevel) {
        ++binsPerLevel[currLevelIdx];
        ++currIdx;
      } else
        break;
    }
    if (currLevelIdx < numLevels - 1)  // boundary case because of numerical pbm
      ++currLevelIdx;
    currLevel += rangeLevel;
  }

  for (int i = 0; i < numLevels; ++i) {
    double minLevel = i * 100 / numLevels;
    double maxLevel = (i + 1) * 100 / numLevels;

    double percentBinsInLevel = 100 * binsPerLevel[i] / double(numWSBins);
    cout << "% bins with " << minLevel << "% - " << maxLevel
         << "% whitespace : " << percentBinsInLevel << endl;
  }
  cout << endl;
  cout << "Median Whitespace % : " << 100 * medianWS << endl;
  cout << "Average Whitespace % : " << 100 * avgWS << endl;
  cout << "Max Whitespace % : " << 100 * maxWS << endl;
  cout << "Min Whitespace % : " << 100 * minWS << endl;

  // now plot in gnuplot format
  gpFile << "#Use this file as a script for gnuplot" << endl;
  gpFile << "#(See http://www.gnuplot.info/ for details)" << endl;
  gpFile << "#" << TimeStamp() << User() << Platform() << endl << endl;
  gpFile << "set title \" " << baseFileName << ", #Cells= " << getNumCells()
         << ", #Bins= " << numWSBins << " \" font \"Times , 22\"" << endl
         << endl;

  gpFile << "set xlabel \" % Whitespace in Bins\" " << endl;
  gpFile << "set ylabel \" % #Bins \" " << endl;
  gpFile << "set xtics 10 " << endl;
  gpFile << "set ytics 10 " << endl;
  gpFile << "set nokey" << endl;

  gpFile << "#   Uncomment these two lines starting with \"set\"" << endl;
  gpFile << "#   to save an EPS file for inclusion into a latex document"
         << endl;
  gpFile << "# set terminal postscript eps color solid 20" << endl;
  gpFile << "# set output \"" << baseFileName << ".eps\"" << endl << endl;

  gpFile << "#   Uncomment these two lines starting with \"set\"" << endl;
  gpFile << "#   to save a PS file for printing" << endl;
  gpFile << "# set terminal postscript portrait color solid 20" << endl;
  gpFile << "# set output \"" << baseFileName << ".ps\"" << endl
         << endl
         << endl;

  gpFile << " plot[0:100][0:100] '-' w l " << endl;
  gpFile << endl << endl;

  for (int i = 0; i < numLevels; ++i) {
    double minLevel = i * 100 / numLevels;
    double maxLevel = (i + 1) * 100 / numLevels;

    double percentBinsInLevel = 100 * binsPerLevel[i] / double(numWSBins);
    gpFile << "  " << minLevel << "  0 " << endl;
    gpFile << "  " << minLevel << "  " << percentBinsInLevel << endl;
    gpFile << "  " << maxLevel << "  " << percentBinsInLevel << endl;
    gpFile << "  " << maxLevel << "  0 " << endl;
    gpFile << endl;
  }
  gpFile << "EOF" << endl << endl;
  gpFile << "pause -1 'Press any key' " << endl;
  gpFile.close();
}

// BEGIN DEPRICATED OLD PLOTTING UTILITY
// KEPT FOR BACKWARD COMPATABILITY
// Functions below here should not be used in creating new plot commands and
// options newer alternatives for all of them exist. Please see RBPlacePlot.h
// for more info on how to modify this file.

void RBPlacement::saveAsPlot(
    const Plotters::RBPlacePlotter::Parameters& plotterParams,
    const char* baseFileName) const {
  std::vector<unsigned> nodes;
  for (unsigned i = 0; i < _hgWDims->getNumNodes(); i++) {
    if (_placement[i].x < plotterParams.xMax &&
        _placement[i].x > plotterParams.xMin &&
        _placement[i].y < plotterParams.yMax &&
        _placement[i].y > plotterParams.yMin)
      nodes.push_back(i);
  }
  plot(plotterParams, nodes, baseFileName);
  nodes.clear();
}

void RBPlacement::plotDummy(
    const Plotters::RBPlacePlotter::Parameters& plotterParams,
    const char* baseFileName) const {
  std::vector<unsigned> nodes;
  for (unsigned i = 0; i < _hgWDims->getNumNodes(); i++) {
    const std::string nodeName = _hgWDims->getNodeNameByIndex(i).c_str();
    if (nodeName.find("dummy") != std::string::npos) nodes.push_back(i);
  }
  plot(plotterParams, nodes, baseFileName);
}

void RBPlacement::plotGrpConstr(
    const Plotters::RBPlacePlotter::Parameters& plotterParams,
    const char* baseFileName) const {
  std::vector<unsigned> nodes;
  for (unsigned i = 0; i < groupRegionConstr.size(); ++i) {
    for (unsigned j = 0; j < groupRegionConstr[i].second.size(); ++j) {
      nodes.push_back(groupRegionConstr[i].second[j]);
    }
  }
  plot(plotterParams, nodes, baseFileName);
}

void RBPlacement::plot(
    const Plotters::RBPlacePlotter::Parameters& plotterParams,
    const std::vector<unsigned>& nodes, const char* baseFileName) const {
  cout << "Saving " << baseFileName << ".gpl" << endl;

  Plotters::RBPlacePlotter plot(plotterParams, baseFileName, *this);
  const Plotters::nodeSet plotNodes(nodes.begin(), nodes.end());

  if (plotterParams.plotNodes) plot.addNodeSet(plotNodes);
  if (plotterParams.plotNodesNames) plot.addNodeSetForLabels(plotNodes);
  if (plotterParams.plotNets) plot.addNodeSetForNets(plotNodes);

  ofstream gpFile;
  gpFile.open((string(baseFileName) + string(".gpl")).c_str());
  gpFile << plot;
  gpFile.close();
}
