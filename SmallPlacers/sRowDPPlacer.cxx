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

// Base Small Placer
//  - Small Placement Problem
//    created by Andrew Caldwell  981206

#include "sRowDPPlacer.h"

#include "stripNetStacks.h"

using std::vector;

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

struct table_compare
// : public binary_function<DPOneRowTableEntry*,DPOneRowTableEntry*,bool>
{
  bool operator()(DPOneRowTableEntry* x, DPOneRowTableEntry* y) {
    return x->_stacks->getCost() > y->_stacks->getCost();
  }
};

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

DPSmallPlacer::DPSmallPlacer(SmallPlacementProblem& problem, Parameters params,
                             MaxMem* maxMem)
    : BaseSmallPlacer(problem, params, maxMem),
      _curSoln(problem.getNumCells()),
      _bestSolnOverPass(problem.getNumCells()) {
  abkfatal(problem.getNumRows() == 1,
           "DP SmallPlacer can't deal with multiple rows");
  abkfatal(problem.getNumCells() >= 2,
           "DP SmallPlacer requires at least 2 cells");

  unsigned bits_in_uns = sizeof(unsigned) * 8;
  abkfatal(problem.getNumNets() < 2 * bits_in_uns,
           "DP SmallPlacer can't deal with this many nets");

  Placement& soln = _problem.getSolnBuffer()[0].placement;
  //  const unsigned numNets  = _problem.getNumNets() ;
  const unsigned numCells = _problem.getNumCells();

  buildConnectivity();

  // set y-coordinates (single row, so these don't change!)
  unsigned n;
  for (n = numCells; n > 0;) {
    --n;
    _curSoln[n].y = soln[n].y = problem.getRows()[0].yMin;
  }

  // for the moment, split the problem into two equal parts.

  bool usingInitSoln = false;

  _q.reserve(numCells);  // workspace to get init soln.

  if (_problem.getSolnBuffer()[0].populated && _params.useInitSoln) {
    for (n = 0; n < numCells; n++) _q.push_back(n);

    std::sort(_q.begin(), _q.end(), DPCompareNodesByInitPl(soln));
    usingInitSoln = true;
  } else {
    if (_params.startSoln == Parameters::SmartStart) {
      initialSoln();  // nb: sorting of _q done in initialSoln()!!
    } else {
      // randomly shuffle the cells...
      for (n = 0; n < numCells; n++) _q.push_back(n);

      RandomUnsigned ru(0, _q.size());
      for (n = 100 * _q.size(); n != 0; n--) {
        unsigned u = ru;
        unsigned v = ru;

        // debug... should rand upper lim be #, or #-1?
        if (u >= numCells) u = numCells - 1;
        if (v >= numCells) v = numCells - 1;

        unsigned tmp = _q[u];
        _q[u] = _q[v];
        _q[v] = tmp;
      }
    }
  }

  // using _q workspace, build two sequences to interleave.
  // list composition depends on the splitType.

  vector<unsigned> _a, _b;
  if (_params.splitType == Parameters::RandomSplit) {
    RandomUnsigned ru(0, _q.size() - 1);

    // continue until no improvement.

    double bestCost = DBL_MAX;  //, currCost = 0. ;
    for (unsigned pass = 0; pass != 20; pass++) {
      _a.clear();
      _b.clear();
      unsigned split = ru;
      if (split >= numCells - 1) split = numCells - 2;

      for (n = 0; n < numCells; n++) {
        if (n <= split)
          _a.push_back(_q[n]);
        else
          _b.push_back(_q[n]);
      }

      double currCost = go(_a, _b);  // do it!
      if (currCost < bestCost - 1.0e-8) {
        bestCost = currCost;
        _problem.getSolnBuffer()[0].populated = true;
        for (n = 0; n != numCells; n++) {
          soln[n].x = _curSoln[n].x;
        }
      } else
        break;
    }
    goto end;
    return;
  } else if (_params.splitType == Parameters::EvenSplit) {
    double bestCost = DBL_MAX;  //, currCost = 0. ;
    for (unsigned pass = 0; pass != 20; pass++) {
      _a.clear();
      _b.clear();
      unsigned split = numCells / 2 - 1;
      for (n = 0; n < numCells; n++) {
        if (n <= split)
          _a.push_back(_q[n]);
        else
          _b.push_back(_q[n]);
      }

      double currCost = go(_a, _b);  // do it!
      if (currCost < bestCost - 1.0e-8) {
        bestCost = currCost;
        _problem.getSolnBuffer()[0].populated = true;
        for (n = 0; n != numCells; n++) {
          soln[n].x = _curSoln[n].x;
        }
      } else
        break;
    }
    goto end;
    return;
  } else if (_params.splitType == Parameters::AllSplit) {
    vector<unsigned> bestOrd(numCells);
    for (n = 0; n != numCells; n++) bestOrd[n] = _q[n];

    double bestCost = DBL_MAX;
    for (unsigned pass = 0; pass != 20; pass++) {
      // consider all splits of the order.
      double bestCostOverPass = DBL_MAX;
      for (unsigned split = 0; split != numCells; split++) {
        _a.clear();
        _b.clear();
        for (n = 0; n != numCells; n++) {
          if (n <= split)
            _a.push_back(bestOrd[n]);
          else
            _b.push_back(bestOrd[n]);
        }
        double currCost = go(_a, _b);

        if (currCost < bestCostOverPass - 1.0e-8) {
          // nb: best soln during this pass; use to prep next.
          bestCostOverPass = currCost;
          for (n = 0; n != numCells; n++) {
            bestOrd[n] = _q[n];
            _bestSolnOverPass[n].x = _curSoln[n].x;
          }
        }
      }
      if (bestCostOverPass < bestCost - 1.0e-8) {
        bestCost = bestCostOverPass;
        _problem.getSolnBuffer()[0].populated = true;
        for (n = 0; n != numCells; n++) soln[n].x = _bestSolnOverPass[n].x;
      } else
        break;
    }
  } else {
    unsigned split = numCells / 2 - 1;  // default to even splitting.
    for (n = 0; n < numCells; n++) {
      if (n <= split)
        _a.push_back(_q[n]);
      else
        _b.push_back(_q[n]);
    }

    go(_a, _b);  // do it!
    goto end;
    return;
  }

end:
  // problem.plot("plot");
  return;
}

////////////////////////////////////////////////////////////////////////////////

DPSmallPlacer::~DPSmallPlacer() {
  if (_adjNets) delete[] _adjNets;
  if (_firstNets) delete[] _firstNets;

  unsigned i;
  for (i = 0; i != _tablePool.size(); i++) {
    delete _tablePool[i]->_stacks;
    delete _tablePool[i];
  }
}

////////////////////////////////////////////////////////////////////////////////

DPOneRowTableEntry* DPSmallPlacer::getTableEntryFromPool() {
  DPOneRowTableEntry* n;
  if (_tablePool.size() != 0) {
    n = _tablePool.back();
    _tablePool.pop_back();
    return n;
  }
  n = new DPOneRowTableEntry;
  n->x[0] = 0;
  n->x[1] = 0;
  n->_stacks = new StripNetStacks(_problem, _firstNets);
  n->parent = 0;
  return n;
}

////////////////////////////////////////////////////////////////////////////////

void DPSmallPlacer::releaseTableEntryToPool(DPOneRowTableEntry* p) {
  p->parent = 0;  // nb: we don't clear anything!!!
  _tablePool.push_back(p);
}

////////////////////////////////////////////////////////////////////////////////

double DPSmallPlacer::go(vector<unsigned>& _a, vector<unsigned>& _b) {
  unsigned i, j;

  // call alternative routing if -ve cells.

  for (i = 0; i != _q.size(); i++) {
    if (_problem.getCellWidth(_q[i]) < 0.0) {
      return go1(_a, _b);
    }
  }

  // routine which does the actual DP.
  const unsigned numCells = _problem.getNumCells();

  //++++ build the (full) DP table.

  unsigned tableDim1 = _a.size() + 1;
  unsigned tableDim2 = _b.size() + 1;

  DPOneRowTableEntry*** table = new DPOneRowTableEntry**[tableDim1];
  for (i = tableDim1; i--;) {
    table[i] = new DPOneRowTableEntry*[tableDim2];

    for (j = tableDim2; j--;) {
      table[i][j] = getTableEntryFromPool();
      table[i][j]->x[0] = i;
      table[i][j]->x[1] = j;
      table[i][j]->_dir = -1;
    }
  }

  //++++ fill in the table; go down the first row and column, then fill
  //     in the rest of the table.

  StripNetStacks worki(_problem, _firstNets);
  StripNetStacks workj(_problem, _firstNets);

  // nb: important!!! we need to "clean" entry (0,0).

  table[0][0]->_stacks->copy(worki);  // worki is clean here.
  table[0][0]->_dir = -1;

  for (i = 1; i < tableDim1; i++) {
    table[i][0]->_stacks->copy(*(table[i - 1][0]->_stacks));
    table[i][0]->_stacks->addNodeToLeft(_a[i - 1]);
    table[i][0]->_dir = 0;
  }
  for (j = 1; j < tableDim2; j++) {
    table[0][j]->_stacks->copy(*(table[0][j - 1]->_stacks));
    table[0][j]->_stacks->addNodeToLeft(_b[j - 1]);
    table[0][j]->_dir = 1;
  }

  for (i = 1; i < tableDim1; i++) {
    unsigned ai = _a[i - 1];

    for (j = 1; j < tableDim2; j++) {
      unsigned bj = _b[j - 1];

      worki.copy(*(table[i - 1][j]->_stacks));
      worki.addNodeToLeft(ai);

      workj.copy(*(table[i][j - 1]->_stacks));
      workj.addNodeToLeft(bj);

      if (worki.getCost() <= workj.getCost()) {
        table[i][j]->_dir = 0;
        table[i][j]->_stacks->copy(worki);
      } else {
        table[i][j]->_dir = 1;
        table[i][j]->_stacks->copy(workj);
      }
    }
  }

  // trace back solution...

  vector<unsigned> order;
  order.reserve(numCells);

  unsigned idx = numCells;

  _maxMem->update("Dynamic Programming Single Row SmallPlacer");

  i = _a.size();
  j = _b.size();

  int dir = table[i][j]->_dir;
  while (dir >= 0) {
    --idx;

    _q[idx] = order[idx] = (dir == 0) ? _a[i - 1] : _b[j - 1];

    if (dir == 0) {
      i--;
    } else if (dir == 1) {
      j--;
    } else
      ;  // huh?

    dir = table[i][j]->_dir;
  }

  worki.copy(*(table[0][0]->_stacks));  // entry (0,0) is always clean.

  _curWL = 0.0;
  for (idx = 0; idx < numCells; idx++) {
    unsigned c = order[idx];

    _curSoln[c].x = worki.getLeftEdge();
    _curWL = worki.addNodeToLeft(c);
  }

  //++++ release allocated memory to the pool.

  for (i = tableDim1; i--;) {
    for (j = tableDim2; j--;) {
      releaseTableEntryToPool(table[i][j]);
    }
    delete[] table[i];
  }
  delete[] table;

  return _curWL;
}

////////////////////////////////////////////////////////////////////////////////

void DPSmallPlacer::buildConnectivity() {
  const SmallPlacementNetlist& netlist = _problem.getNetlist();

  unsigned numCells = _problem.getNumCells();
  unsigned numNets = _problem.getNumNets();
  unsigned numPins = numNets * (numCells + 1);

  _adjNets = new EdgeAndOffset[numPins];
  _firstNets = new EdgeAndOffset*[numCells + 1];

  unsigned c, n, idx = 0;

  // nets (and offsets) adjacent to each node.
  for (c = 0; c < _problem.getNumCells(); c++) {
    _firstNets[c] = _adjNets + idx;
    for (n = 0; n < _problem.getNumNets(); n++) {
      const BBox& pinOffset = netlist.getPinOffset(c, n);
      if (!pinOffset.isEmpty()) {
        _adjNets[idx].edge = n;
        _adjNets[idx].pinOffset = pinOffset;
        idx++;
      }
    }
  }
  _firstNets[c] = _adjNets + idx;
}

////////////////////////////////////////////////////////////////////////////////

void DPSmallPlacer::initialSoln() {
  const SmallPlacementNetlist& netlist = _problem.getNetlist();

  unsigned numCells = _problem.getNumCells();
  unsigned numNets = _problem.getNumNets();
  unsigned numPins = numNets * (numCells + 1);

  double xMin = _problem.getRows()[0].xMin;
  double xMax = _problem.getRows()[0].xMax;

  unsigned* edgeDeg = new unsigned[numNets];
  unsigned* cellDeg = new unsigned[numCells];

  double* fixedTerminal = new double[numNets];

  unsigned* adjCell = new unsigned[numPins + 1];
  unsigned** firstCell = new unsigned*[numNets + 1];

  unsigned* adjNet = new unsigned[numPins + 1];
  unsigned** firstNet = new unsigned*[numCells + 1];

  unsigned n, c, idx;

  // determine edges adjacent to cells.
  for (idx = 0, n = 0; n < numNets; n++) {
    edgeDeg[n] = 0;

    firstCell[n] = adjCell + idx;
    for (c = 0; c < numCells; c++) {
      // double pinOffset = netlist.getPinOffset(c,n).x ;
      // if (pinOffset != DBL_MAX)
      if (!netlist.getPinOffset(c, n).isEmpty()) {
        adjCell[idx++] = c;
        edgeDeg[n]++;
      }
    }

    const BBox& netBBox = _problem.getNetTerminalBBox(n);
    if (!netBBox.isEmpty()) {
      fixedTerminal[n] = (netBBox.xMin + netBBox.xMax) / 2.;
      edgeDeg[n]++;
    } else
      fixedTerminal[n] = 0.;
  }
  firstCell[n] = adjCell + idx;

  for (idx = 0, c = 0; c < numCells; c++) {
    cellDeg[c] = 0;

    firstNet[c] = adjNet + idx;
    for (n = 0; n < numNets; n++) {
      // double pinOffset = netlist.getPinOffset(c,n).x ;
      // if (pinOffset != DBL_MAX)
      if (!netlist.getPinOffset(c, n).isEmpty()) {
        adjNet[idx++] = n;
        cellDeg[c]++;
      }
    }

    _curSoln[c].x = (xMin + xMax) / 2.;
  }
  firstNet[c] = adjNet + idx;

  // perform cog iterations.
  for (unsigned iter = 1; iter != 11; iter++) {
    for (unsigned curr = 0; curr < numCells; curr++)  // update cell i posn.
    {
      double cellPosn = 0.;
      for (unsigned* i = firstNet[curr]; i != firstNet[curr + 1]; i++) {
        n = *i;

        double netPosn = fixedTerminal[n];
        for (unsigned* j = firstCell[n]; j != firstCell[n + 1]; j++) {
          c = *j;
          // netPosn += (_curSoln[c].x + netlist.getPinOffset(c,n).x) ;
          const BBox& pinOffset = netlist.getPinOffset(c, n);
          netPosn += (_curSoln[c].x + 0.5 * (pinOffset.xMin + pinOffset.xMax));
        }
        cellPosn += (netPosn / ((double)edgeDeg[n]));
      }
      if (cellDeg[curr] != 0)
        _curSoln[curr].x = cellPosn / ((double)cellDeg[curr]);
    }
  }

  for (c = 0; c < numCells; c++) _q.push_back(c);
  std::sort(_q.begin(), _q.end(), DPCompareNodesByInitPl(_curSoln));

  delete[] adjNet;
  delete[] firstNet;
  delete[] adjCell;
  delete[] firstCell;

  delete[] fixedTerminal;
  delete[] edgeDeg;
  delete[] cellDeg;
}

////////////////////////////////////////////////////////////////////////////////

double DPSmallPlacer::go1(vector<unsigned>& _a, vector<unsigned>& _b) {
  const unsigned numCells = _problem.getNumCells();

  //++++++++++++ new implementation starts here.

  double left_edge = _problem.getRows()[0].xMin;

  // strip out negative WS cells.

  unsigned i;
  vector<unsigned> cells[3];  // lists of cells; third is negWS cells.

  for (i = 0; i != _a.size(); i++) {
    double width = _problem.getCellWidth(_a[i]);
    if (width >= 0.) {
      cells[0].push_back(_a[i]);
    } else {
      cells[2].push_back(_a[i]);
    }
  }
  for (i = 0; i != _b.size(); i++) {
    double width = _problem.getCellWidth(_b[i]);
    if (width >= 0.) {
      cells[1].push_back(_b[i]);
    } else {
      cells[2].push_back(_b[i]);
    }
  }

  // establish dimensions and offsets for keys.

  unsigned dim[3], off[3];

  dim[0] = cells[0].size() + 1;
  dim[1] = cells[1].size() + 1;
  dim[2] = cells[2].size() + 1;

  off[0] = 0;
  off[1] = dim[0];
  off[2] = dim[1] * off[1];

  unsigned y[3];
  unsigned max_index = dim[0] * dim[1] * dim[2];

  // keep track of what's in the open queue.

  vector<unsigned> visited;  // keys (states) which are visited during the
                             // search. used for releasing allocated memory.
  visited.reserve(max_index);  // never more than this.

  vector<bool> in_open(max_index, false);
  vector<bool> in_clos(max_index,
                       false);  // don't need to maintain closed queue.
  vector<DPOneRowTableEntry*> open_queue;  // open queue; impl. as a heap.

  StripNetStacks init(_problem, _firstNets);  // clean copy of stacks.
  StripNetStacks work(_problem, _firstNets);

  // simple for now... just allocate every possible entry.
  vector<DPOneRowTableEntry*> all_entries(max_index, 0);

  // get start state (which has key = 0).

  all_entries[0] = getTableEntryFromPool();
  visited.push_back(0);

  all_entries[0]->x[0] = 0;
  all_entries[0]->x[1] = 0;
  all_entries[0]->x[2] = 0;
  all_entries[0]->parent = 0;
  all_entries[0]->_stacks->copy(init);  // nb: make sure we're clean.

  // prep the open queue.
  open_queue.push_back(all_entries[0]);
  in_open[0] = true;
  std::make_heap(open_queue.begin(), open_queue.end(), table_compare());

  // perform the search.

  DPOneRowTableEntry* best_table_entry = 0;
  while (open_queue.size()) {
    // get cheapest element.
    DPOneRowTableEntry* first = open_queue.front();
    std::pop_heap(open_queue.begin(), open_queue.end(), table_compare());
    open_queue.pop_back();

    unsigned* x = first->x;

    unsigned key = x[0] + x[1] * off[1] + x[2] * off[2];

    in_open[key] = false;
    in_clos[key] = true;

    // finished?
    unsigned num_placed = first->x[0] + first->x[1] + first->x[2];

    if (num_placed == numCells) {
      double curr_cost = first->_stacks->getCost();
      best_table_entry = all_entries[key];
      break;
    } else {
      // generate children; from each list into row.

      for (unsigned lst = 0; lst != 3; lst++) {
        unsigned next = first->x[lst];

        unsigned leng = cells[lst].size();

        if (next >= leng) continue;  // current list exhausted.

        unsigned cell = cells[lst][next];

        // compute cost of adding cell.

        work.copy(*(first->_stacks));  // use the 'work' stack made at
                                       // beginning.
        work.addNodeToLeft(cell);

        y[0] = x[0];
        y[1] = x[1];
        y[2] = x[2];
        y[lst] += 1;

        key = y[0] + y[1] * off[1] + y[2] * off[2];

        if (all_entries[key] == 0) {
          all_entries[key] = getTableEntryFromPool();
          abkassert(all_entries[key] != 0 && all_entries[key]->_stacks != 0,
                    " error: ask Andrew Kennings ");
          visited.push_back(key);

          all_entries[key]->x[0] = y[0];
          all_entries[key]->x[1] = y[1];
          all_entries[key]->x[2] = y[2];
          all_entries[key]->_stacks->copy(
              init);  // nb: make sure we're clean. need this here???
        } else {
          abkfatal(all_entries[key] != 0, "internal problem; 1 row DP");
        }

        DPOneRowTableEntry* child = all_entries[key];

        // following bounds away from negative cells going first.
        if (work.getLeftEdge() < left_edge) {
          continue;
        }

        // see if insersion into the queue is required.

        if (!in_open[key] && !in_clos[key]) {
          all_entries[key]->parent = first;
          all_entries[key]->_stacks->copy(work);

          // cout<<"------>new insersion; cost: "
          //     <<all_entries[key]->_stacks->getCost()<<endl;

          open_queue.push_back(all_entries[key]);
          std::push_heap(open_queue.begin(), open_queue.end(), table_compare());
          in_open[key] = true;
        } else {
          if (in_open[key]) {
            // cout<<"------>currently in open..."<<endl;

            // nb: if in the open queue, then the cost stored in the keyed
            //     entry will be it's cost in the heap.

            double newCost = work.getCost();
            double curCost = all_entries[key]->_stacks->getCost();

            if (newCost < curCost || (newCost <= curCost && lst == 0)) {
              // find in the heap.
              // cout<<"------>replace... new cost is better:
              // "<<all_entries[key]->_stacks->getCost()<<endl;

              vector<DPOneRowTableEntry*>::iterator loc;
              for (loc = open_queue.begin(); loc != open_queue.end(); loc++) {
                if ((*loc) == all_entries[key]) {
                  break;
                }
              }

              (*loc)->parent = first;
              (*loc)->_stacks->copy(work);
              std::make_heap(open_queue.begin(), open_queue.end(),
                             table_compare());
            } else {
              // cout<<"------>ignore... current cost
              // ("<<all_entries[key]->_stacks->getCost()<<")"
              //     <<" is better than new cost: "<<work.getCost()<<endl;
              ;
            }
          } else {
            // cout<<"------>currently in closed..."<<endl;

            // currently in closed queue; means children have been expanded...
            // re-insersion of this node into open (if cheaper) will be revisit
            // all children too.

            double newCost = work.getCost();
            double curCost = all_entries[key]->_stacks->getCost();
            if (newCost < curCost) {
              all_entries[key]->parent = first;
              all_entries[key]->_stacks->copy(work);

              // cout<<"------>re-insersion: cost:
              // "<<all_entries[key]->_stacks->getCost()<<endl;

              open_queue.push_back(all_entries[key]);
              std::push_heap(open_queue.begin(), open_queue.end(),
                             table_compare());
              in_open[key] = true;
              in_clos[key] = false;
            } else {
              // cout<<"------>ignore... previous cost better:
              // "<<all_entries[key]->_stacks->getCost()<<endl;
              ;
            }
          }
        }
      }
    }
  }

  // trace back solution here....

  double retval = DBL_MAX;
  if (best_table_entry) {
    abkassert(best_table_entry->x[0] + best_table_entry->x[1] +
                      best_table_entry->x[2] ==
                  numCells,
              " error: ask Andrew Kennings ");

    vector<DPOneRowTableEntry*> traceback;

    DPOneRowTableEntry* curr = best_table_entry;
    while (curr) {
      traceback.push_back(curr);
      curr = curr->parent;
    }

    curr = traceback.back();

    unsigned idx[3];
    idx[0] = 0;
    idx[1] = 0;
    idx[2] = 0;

    StripNetStacks work(_problem, _firstNets);

    unsigned cnt = 0;
    for (i = traceback.size(); i--;) {
      curr = traceback[i];

      if (i != 0) {
        // determine which coordinate differs.

        unsigned lst = 0;

        unsigned coord;
        for (coord = 0; coord != 3; coord++) {
          if (traceback[i]->x[coord] != traceback[i - 1]->x[coord]) {
            lst = coord;
            break;
          }
        }

        unsigned cell = cells[lst][idx[lst]++];

        _q[cnt++] = cell;
        _curSoln[cell].x = work.getLeftEdge();

        retval = work.addNodeToLeft(cell);
      }
    }
  }

  // cleanup.

  for (i = 0; i < visited.size(); i++) {
    releaseTableEntryToPool(all_entries[visited[i]]);
  }

  return retval;
}

////////////////////////////////////////////////////////////////////////////////
