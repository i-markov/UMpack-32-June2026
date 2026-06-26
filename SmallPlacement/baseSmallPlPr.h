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

// Small Placement Interface
//  - Small Placement Problem
//    created by Andrew Caldwell  981206
//  CHANGES
//   990501 ilm renamed into BaseSmallPlacementProblem
//   001128 s.a added a plotter and a fnuction to detect/remove overlaps

#ifndef __BASE_SMALLPL_PROB_H__
#define __BASE_SMALLPL_PROB_H__

#include <algorithm>
#include <cfloat>
#include <fstream>
#include <vector>

#include "ABKCommon/abkassert.h"
#include "ABKCommon/abktypes.h"
#include "ABKCommon/verbosity.h"
#include "Ctainers/listO1size.h"
#include "Placement/symOri.h"
#include "smallNetlist.h"
#include "smallPlRow.h"

class HGraphWDimensions;

struct SmallPlacementSolution {
  struct CompCellsByLoc {
    const Placement& pl;
    CompCellsByLoc(const Placement& place) : pl(place) {}
    bool operator()(unsigned n1, unsigned n2) const {
      return (pl[n1].x < pl[n2].x);
    }
  };

  Placement placement;
  std::vector<Orient> orientations;
  bool populated;

  SmallPlacementSolution(unsigned numCells = 0)
      : placement(numCells), orientations(numCells), populated(false) {}

  // sets perm to be a permutation of 1..n (where n == num cells)
  // coresponding to their left->right ordering in the placement
  void setPermutation(std::vector<unsigned>& perm) const {
    perm.clear();
    perm.reserve(placement.getSize());
    for (unsigned n = 0; n < placement.getSize(); n++) perm.push_back(n);
    std::sort(perm.begin(), perm.end(), CompCellsByLoc(placement));
  }
};

class BaseSmallPlacementProblem {
 protected:
  SmallPlacementNetlist _netlist;
  // gives the connections, and pin-offsets, for each cell-net pair

  std::vector<SmallPlacementRow> _rows;
  // note: row contains site information. However,
  // first generation(ISPD) placers will not use it.
  // They will assume that:
  //  1) the entier width of the row is available for cells
  //  2) white space (even negative ws) is to be equally
  //	distributed

  bit_vector _fixed;
  // this gives the fixed cells in instance. by sadya

  std::vector<double> _cellWidths;
  // all cells are assumed to be exactly 1 row high

  std::vector<BBox> _netTerminalBBoxes;
  // for each net, this gives the BBox of its terminal's pins

  std::vector<SmallPlacementSolution> _solnBuffer;
  unsigned _bestSoln;
  double _cost;

  Verbosity _verb;

  void readAUX(const char* auxFile);
  bool readPLN(std::ifstream& plnStream);
  bool readPLC(std::ifstream& plcStream);
  bool readPLA(std::ifstream& plcStream);

  bool savePLN(std::ofstream& plnStream);
  bool savePLC(std::ofstream& plcStream);
  bool savePLA(std::ofstream& plaStream);

  // used by derived classes
  BaseSmallPlacementProblem(unsigned numCells, Verbosity verb)
      : _bestSoln(0), _cost(-DBL_MAX), _verb(verb) {
    _solnBuffer.push_back(SmallPlacementSolution(numCells));
  }

  // used by derived classes
  BaseSmallPlacementProblem(unsigned numCells, unsigned numNets, Verbosity verb)
      : _netlist(numCells, numNets),
        _bestSoln(0),
        _cost(-DBL_MAX),
        _verb(verb) {
    _solnBuffer.push_back(SmallPlacementSolution(numCells));
  }

 public:
  BaseSmallPlacementProblem(const char* auxFile, Verbosity verb = Verbosity())
      : _bestSoln(0), _cost(-DBL_MAX), _verb(verb) {
    readAUX(auxFile);
    _solnBuffer.push_back(SmallPlacementSolution(_netlist.getNumCells()));
  }

  virtual ~BaseSmallPlacementProblem() {}

  const SmallPlacementNetlist& getNetlist() const { return _netlist; }
  const std::vector<double>& getCellWidths() const { return _cellWidths; }
  double getCellWidth(unsigned c) const {
    abkassert(c < _netlist.getNumCells(), "out of bounds c");
    return _cellWidths[c];
  }

  void save(const char* baseFileName);

  std::vector<SmallPlacementSolution>& getSolnBuffer() { return _solnBuffer; }
  SmallPlacementSolution& getBestSoln() { return _solnBuffer[_bestSoln]; }
  unsigned getBestSolnNum() { return _bestSoln; }

  unsigned getNumRows() const { return _rows.size(); }
  const std::vector<SmallPlacementRow>& getRows() const { return _rows; }

  unsigned getNumCells() const { return _netlist.getNumCells(); }
  unsigned getNumNets() const { return _netlist.getNumNets(); }

  const BBox& getNetTerminalBBox(unsigned n) const {
    return _netTerminalBBoxes[n];
  }
  double getCost() const {
    abkfatal(_cost != -DBL_MAX, " Uninitialized cost");
    return _cost;
  }
  void setCost(double cost) { _cost = cost; }

  friend std::ostream& operator<<(std::ostream& os,
                                  const BaseSmallPlacementProblem& prob);
  // added by sadya to plot small placement solutions
  void plot(const char* baseFileName);

  // removes overlaps in the current placement instance if any
  // this function should not be used
  // instead use a similar function in RBPlace to remove overlaps
  void remOverlaps(void);

  // this function returns the wirelength of the placement
  double calculateWL(const Placement& placement);

  void calculateMinXPinOffsets(void) { _netlist.calculateMinXPinOffsets(); }
};
#endif
