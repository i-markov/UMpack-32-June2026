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

#ifndef ROW_BASED_PLACEMENT_PLOT_H
#define ROW_BASED_PLACEMENT_PLOT_H

#include <iostream>
#include <vector>

#include <Geoms/bbox.h>

// forward declaration
class RBPlacement;

namespace Plotters {
typedef unsigned nodeIdxT;
typedef std::vector<nodeIdxT> nodeSet;

enum Colors {
  blackdotted = 0,
  red = 1,
  green = 2,
  blue = 3,
  magenta = 4,
  cyan = 5,
  sienna = 6,
  orange = 7,
  coral = 8
};

// HOW TO PROPERLY MODIFY
//   The RBPlacePlotter class
//	If it is your intention to add a new plot style parameter or otherwise:
//   1.  Create a parameter in the RBPlacePlotter::Parameters struct/class
//   2.  Add default initialization to your parameter in the default constructor
//   "RBPlacePlotter::Parameters::Parameters(void)"
//			Most likely you want it off, please test rigorously if
//you want it on by default
//   3.  Add a command line parameter initialization if desired in the
//   constructor "RBPlacePlotter::Parameters::Parameters(int argc, char **
//   argv)"
//   4.  Create a member function called writeGnuplotMyParameter() or something
//   similar in the RBPlacePlotter class.
//			This function should take as a first parameter "ostream&
//os".
//   5.  Specify the number of data files this function will write getFileCount
//   function (this may be 0).
//   6.  Add a section to writeDataFileContents() containing a call to your
//   function (writeGnuplotMyParameter())
//	7.  If your parameter requires changing the file count, modify the
//getFileCount() function, and
//           be sure to call call the update() function after your data is set
//           (for an example see the addNodeSet() function)
//   8.  If your parameter is difficult to use, consider adding a static member
//   function that knows how to set up and initialize it.
//			See RBPlacePlotter::plotFreeShapeFloorPlan() as an
//example.
//   9.  You should not have to modify the operator<< for RBPlacePlotter, but
//   please look there to make sure your
//			Parameter plays nicely with this function
//   10.  As a backup check, Search for "plotNodes" in RBPlacePlot.cxx and
//   RBPlacePlot.h
//			be sure you have mimicked all code for your new
//parameter and function

class RBPlacePlotter {
  friend std::ostream& operator<<(std::ostream& os, const RBPlacePlotter& rbp);

 public:
  struct Parameters {
    Parameters(void);
    Parameters(const Parameters&);
    Parameters(int argc, const char** argv, bool = true);

    bool plotNodes;
    bool plotNodesNames;
    bool plotNets;
    bool plotSteinerNets;
    bool plotSiteMap;
    bool plotRows;
    bool plotFPOutlines;
    bool plotCutLines;
    bool plotFillerCells;
    bool oneColorPerGroup;

    double xMin;
    double xMax;
    double yMin;
    double yMax;

    std::string FPOutlinesName;
    std::string CutLinesName;
  };

  RBPlacePlotter(RBPlacePlotter::Parameters modes, std::string baseFileName,
                 const RBPlacement& rbp);

  void addNodeSet(const nodeSet& ns);
  void addNodeSetForLabels(const nodeSet& ns);
  void addNodeSetForNets(const nodeSet& ns);

  void assignRandomColors(void) const;
  void assignStandardColors(void) const;
  static void plotFreeshapeFloorplan(std::string baseFileName,
                                     const RBPlacement& rbp, int argc = 0,
                                     const char** argv = 0);
  static void plotCellGroups(std::string baseFileName, const RBPlacement& rbp,
                             unsigned numCellGroups);

 private:
  // Private Member Functions
  void writeGnuplotFileHeader(std::ostream& os) const;
  void writeGnuplotDataFilesNameLabels(std::ostream& os) const;
  void writeGnuplotPlotCommand(std::ostream& os) const;
  void writeGnuplotDataFilesContents(std::ostream& os) const;
  void writeGnuplotFileFooter(std::ostream& os) const;

  void writeGnuplotLargeFixedNodes(std::ostream& os,
                                   const nodeSet& nodes) const;
  void writeGnuplotLargeMovableNodes(std::ostream& os,
                                     const nodeSet& nodes) const;
  // void writeGnuplotLargeNodes(std::ostream& os, const nodeSet& nodes) const;
  // void writeGnuplotSmallNodes(std::ostream& os, const nodeSet& nodes) const;
  void writeGnuplotSmallFixedNodes(std::ostream& os,
                                   const nodeSet& nodes) const;
  void writeGnuplotSmallMovableNodes(std::ostream& os,
                                     const nodeSet& nodes) const;

  void writeGnuplotNodeNameLabels(std::ostream& os, const nodeSet& nodes) const;
  void writeGnuplotNets(std::ostream& os, const nodeSet& nodes) const;
  void writeGnuplotSteinerNets(std::ostream& os, const nodeSet& nodes) const;
  void writeGnuplotSiteMap(std::ostream& os) const;
  void writeGnuplotRows(std::ostream& os) const;
  void writeNodesSetsRandomColors(std::ostream& os) const;

  // Call update any time the number
  // of "gnuplot files" (i.e. _file_ct) changes
  void update(void);
  void setFilesAndColors(void);

  size_t getFileCount(void) const;

  // Private Data

  RBPlacePlotter::Parameters _modes;
  std::string _baseFileName;
  const RBPlacement& _rbp;
  std::vector<nodeSet> _nodeGroups;
  std::vector<nodeSet> _nodes2Label;
  std::vector<nodeSet> _nodes4Nets;
  mutable std::vector<Colors> _colors;
  BBox _plotBox;
  size_t _file_ct;
};

}  // namespace Plotters

#endif
