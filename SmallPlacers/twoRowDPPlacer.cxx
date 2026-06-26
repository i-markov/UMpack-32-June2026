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

#include "twoRowDPPlacer.h"

#include <algorithm>
#include <list>

using std::cout;
using std::endl;
using std::make_heap;
using std::max;
using std::pop_heap;
using std::push_heap;
using std::vector;

////////////////////////////////////////////////////////////////////////////////

twoRowDPSmallPlacer::twoRowDPSmallPlacer(SmallPlacementProblem& problem,
                                         Parameters params, MaxMem* maxMem)
    : BaseSmallPlacer(problem, params, maxMem),
      _curSoln(problem.getNumCells()) {
  abkfatal(problem.getNumRows() == 2, "DP2D requires two rows");
  abkfatal(problem.getNumCells() >= 2, "DP2D requires >= 2 cells");

  unsigned bits_in_uns = sizeof(unsigned) * 8;
  abkfatal(problem.getNumNets() < 2 * bits_in_uns,
           "DP2D can't deal with this many nets");

  _problem.calculateMinXPinOffsets();

  Placement& soln = _problem.getSolnBuffer()[0].placement;
  const unsigned numCells = _problem.getNumCells();

  buildConnectivity();

  // create cell list either by random shuffle, or using provided soln.

  unsigned n;
  _usingInitSoln = false;

  _q.reserve(numCells);  // workspace to get init soln.

  _minRowViol[0] = DBL_MAX;
  _minRowViol[1] = DBL_MAX;

  vector<vector<int> > _cellsInRow;
  _cellsInRow.resize(2);

  if (_problem.getSolnBuffer()[0].populated && _params.useInitSoln)
  // if(_params.useInitSoln)
  {
    // Push cells into row 0 until full, then row 1.
    for (n = 0; n < numCells; n++) _q.push_back(n);

    double left[2];
    left[0] = _problem.getRows()[0].xMin;
    left[1] = _problem.getRows()[1].xMin;
    unsigned row = 0;
    for (n = 0; n < numCells; n++) {
      double width = _problem.getCellWidth(_q[n]);
      if (row == 0) {
        // nb: make sure at least one cell in row 0!
        if ((_cellsInRow[0].size() == 0) ||
            (left[0] + width <= _problem.getRows()[0].xMax))
          row = 0;
        else
          row = 1;
      }

      left[row] += width;

      _cellsInRow[row].push_back(_q[n]);
    }

    _usingInitSoln = true;
  } else {
    // If no initial solution, then use something random.
    //++++++++++++ can we be smarter here???
    if (_params.startSoln == Parameters::RandomStart) {
      // randomly shuffle the cells...
      for (n = 0; n < numCells; n++) _q.push_back(n);

      RandomUnsigned ru(0, _q.size());
      for (n = 100 * _q.size(); n != 0; n--) {
        unsigned u = ru;
        unsigned v = ru;

        if (u >= numCells) u = numCells - 1;
        if (v >= numCells) v = numCells - 1;

        unsigned tmp = _q[u];
        _q[u] = _q[v];
        _q[v] = tmp;
      }
    } else {
      for (n = 0; n < numCells; n++) _q.push_back(n);
    }

    double left[2];
    left[0] = _problem.getRows()[0].xMin;
    left[1] = _problem.getRows()[1].xMin;
    unsigned row = 0;
    for (n = 0; n < numCells; n++) {
      double width = _problem.getCellWidth(_q[n]);
      if (row == 0) {
        // nb: make sure at least one cell in row 0!
        if ((_cellsInRow[0].size() == 0) ||
            (left[0] + width <= _problem.getRows()[0].xMax))
          row = 0;
        else
          row = 1;
      }

      left[row] += width;

      _cellsInRow[row].push_back(_q[n]);
    }
  }

  // Determine cost of the initial solution.
  twoRowDPStacks initStk(_problem, _firstNets);

  unsigned r;
  for (r = 0; r <= 1; r++) {
    for (n = 0; n < _cellsInRow[r].size(); n++) {
      int cellId = _cellsInRow[r][n];

      soln[cellId].x = _curSoln[cellId].x = initStk.getLeftEdge(r);
      soln[cellId].y = _curSoln[cellId].y = initStk.getBottomEdge(r);

      initStk.addNode(cellId, r);
    }
  }
  _problem.getSolnBuffer()[0].populated = true;

  _bestCost = initStk.getCost();
  _minRowViol[0] =
      max(0.0, (initStk.getRiteEdge(0) - _problem.getRows()[0].xMax)) +
      max(0.0, (_problem.getRows()[0].xMin - initStk.getLeftEdge(0)));
  _minRowViol[1] =
      max(0.0, (initStk.getRiteEdge(1) - _problem.getRows()[1].xMax)) +
      max(0.0, (_problem.getRows()[1].xMin - initStk.getLeftEdge(1)));

  // cout<<"viol[0]: "<<_minRowViol[0]<<" viol[1]: "<<_minRowViol[1]<<endl ;

  // repeat until no improvement.

  unsigned i = 0;
  double newcost = go(_cellsInRow);
  while (newcost != DBL_MAX) {
    twoRowDPStacks stk(_problem, _firstNets);

    for (r = 0; r <= 1; r++) {
      for (n = 0; n < _cellsInRow[r].size(); n++) {
        int cellId = _cellsInRow[r][n];

        soln[cellId].x = _curSoln[cellId].x;
        soln[cellId].y = _curSoln[cellId].y;

        stk.addNode(cellId, r);
      }
    }
    _problem.getSolnBuffer()[0].populated = true;

    _bestCost = stk.getCost();

    abkfatal(fabs(_bestCost - newcost) < 1e-6,
             "best cost from algorithm"
             " and newly calculated value disagree");

    _minRowViol[0] =
        max(0.0, (stk.getRiteEdge(0) - _problem.getRows()[0].xMax)) +
        max(0.0, (_problem.getRows()[0].xMin - stk.getLeftEdge(0)));
    _minRowViol[1] =
        max(0.0, (stk.getRiteEdge(1) - _problem.getRows()[1].xMax)) +
        max(0.0, (_problem.getRows()[1].xMin - stk.getLeftEdge(1)));

    newcost = go(_cellsInRow);

    // cout<<"cost: "<<_bestCost<<" viol[0]: "<<_minRowViol[0]<<" viol[1]:
    // "<<_minRowViol[1]<<endl ;

    // char buf[256] ;
    // sprintf(buf,"gplot_%d", i) ;
    // i++;
    //_problem.plot(buf);
  }

  //_problem.plot("plot");
}

////////////////////////////////////////////////////////////////////////////////

twoRowDPSmallPlacer::~twoRowDPSmallPlacer() {
  if (_adjNets) delete[] _adjNets;
  if (_firstNets) delete[] _firstNets;

  unsigned i;
  for (i = 0; i != _tablePool.size(); i++) {
    delete _tablePool[i]->st;
    delete _tablePool[i];
  }
}

////////////////////////////////////////////////////////////////////////////////

tableEntry* twoRowDPSmallPlacer::getTableEntryFromPool() {
  tableEntry* n;
  if (_tablePool.size() != 0) {
    n = _tablePool.back();
    _tablePool.pop_back();
    return n;
  }
  n = new tableEntry;
  n->x[0] = 0;
  n->x[1] = 0;
  n->x[2] = 0;
  n->x[3] = 0;
  n->x[4] = 0;
  n->x[5] = 0;
  n->x[6] = 0;
  n->x[7] = 0;

  n->st = new twoRowDPStacks(_problem, _firstNets);
  n->parent = 0;
  return n;
}

////////////////////////////////////////////////////////////////////////////////

void twoRowDPSmallPlacer::releaseTableEntryToPool(tableEntry* p) {
  p->parent = 0;  // nb: we don't clear anything!!!
  _tablePool.push_back(p);
}

////////////////////////////////////////////////////////////////////////////////

void twoRowDPSmallPlacer::dump_open_queue(vector<tableEntry*>& op) {
  vector<tableEntry*>::iterator i;
  for (i = op.begin(); i != op.end(); i++) {
    tableEntry* p = *i;
    for (unsigned c = 0; c != 8; c++) cout << p->x[c] << " ";
    cout << "->cost: " << p->st->getCost() << endl;
  }
  cout << "-------------------------------------" << endl;
}

double twoRowDPSmallPlacer::go(vector<vector<int> >& _list) {
  const unsigned numCells = _problem.getNumCells();

  //++++++++++++ new implementation starts here.

  // record the violation of the best solution found over all iterations;
  // we never want to return something more infeasible, or of worst cost
  // and this helps us bound.

  double left_edge[2], rite_edge[2], best_viol[2], best_cost = _bestCost;

  best_viol[0] = _minRowViol[0];
  best_viol[1] = _minRowViol[1];

  left_edge[0] = _problem.getRows()[0].xMin;
  left_edge[1] = _problem.getRows()[1].xMin;

  rite_edge[0] = _problem.getRows()[0].xMax;
  rite_edge[1] = _problem.getRows()[1].xMax;

  // separate positive and negative width cells.

  vector<vector<int> > posWidthCells(2);
  vector<vector<int> > negWidthCells(2);

  unsigned i, r, l;
  for (r = 0; r != 2; r++) {
    for (i = 0; i < _list[r].size(); i++) {
      int cell = _list[r][i];

      if (_problem.getCellWidth(cell) > 0.0) {
        posWidthCells[r].push_back(cell);
      } else {
        negWidthCells[r].push_back(cell);
      }
    }
  }

  // with negative width cells, the search space is 8 dimensional... :-(

  unsigned dim[8], off[8];

  dim[0] = dim[1] = posWidthCells[0].size() + 1;
  dim[2] = dim[3] = posWidthCells[1].size() + 1;
  dim[4] = dim[5] = negWidthCells[0].size() + 1;
  dim[6] = dim[7] = negWidthCells[1].size() + 1;

  off[0] = 0;
  off[1] = dim[0];
  for (i = 2; i != 8; i++) {
    off[i] = off[i - 1] * dim[i - 1];
  }

  unsigned max_index = 1;
  for (i = 0; i != 8; i++) {
    max_index *= dim[i];
  }

  // keep track of what's in the open queue.

  vector<unsigned> visited;  // keys (states) which are visited during the
                             // search. used for releasing allocated memory.
  visited.reserve(max_index);  // never more than this.

  vector<bool> in_open(max_index, false);
  vector<bool> in_clos(max_index,
                       false);     // don't need to maintain closed queue.
  vector<tableEntry*> open_queue;  // open queue; impl. as a heap.

  twoRowDPStacks init(_problem, _firstNets);  // clean copy of stacks.
  twoRowDPStacks work1(_problem, _firstNets);

  // simple for now... just allocate every possible entry.
  vector<tableEntry*> all_entries(max_index, 0);

  // get start state (which has key = 0).

  all_entries[0] = getTableEntryFromPool();
  visited.push_back(0);

  for (i = 0; i != 8; i++) all_entries[0]->x[i] = 0;
  all_entries[0]->parent = 0;
  all_entries[0]->st->copy(init);  // nb: make sure we're clean.

  // prep the open queue.
  open_queue.push_back(all_entries[0]);
  in_open[0] = true;
  make_heap(open_queue.begin(), open_queue.end(), two_dp_table_compare());
  // dump_open_queue(open_queue);

  // perform the search.

  tableEntry* best_table_entry = 0;

  while (open_queue.size()) {
    // get cheapest element.
    tableEntry* first = open_queue.front();
    pop_heap(open_queue.begin(), open_queue.end(), two_dp_table_compare());
    open_queue.pop_back();

    unsigned* x = first->x;

    unsigned key = x[0];
    for (i = 1; i != 8; i++) key += (x[i] * off[i]);

    // abkfatal(in_open[key] == true, "DP2D internal error") ;
    in_open[key] = false;
    in_clos[key] = true;

    // debug... trace back from current to confirm cost.

    // finished?
    unsigned num_placed = 0;
    for (i = 0; i != 8; i++) {
      num_placed += first->x[i];
    }

    if (num_placed == numCells) {
      double curr_cost = first->st->getCost();

      // nb: if the cost is actually same or worse than that passed into the
      // routine,
      //     break without marking this (best) solution since we have failed to
      //     find something better.

      if (curr_cost >= (best_cost - 1.0e-3)) break;  // no better!

      best_table_entry = all_entries[key];
      break;
    } else {
      // generate children; actually 4 lists to consider now due to pos/neg
      // cells.

      for (unsigned lst = 0; lst != 4; lst++) {
        unsigned indx = lst * 2;
        unsigned next =
            first->x[indx] + first->x[indx + 1];  // index of next cell.

        unsigned leng =
            ((lst == 0 || lst == 1) ? posWidthCells[lst - 0].size()
                                    : negWidthCells[lst - 2].size());

        if (next >= leng) {
          continue;
        }  // current list done.

        unsigned cell = ((lst == 0 || lst == 1) ? posWidthCells[lst - 0][next]
                                                : negWidthCells[lst - 2][next]);

        // cout << "next: " <<next << " leng: " << leng << " cell: " << cell <<
        // endl ;

        for (unsigned row = 0; row != 2; row++) {
          // compute the appropriate child coordinate.

          // compute cost of adding cell.

          work1.copy(*(first->st));  // use the 'work' stack made at beginning.
          work1.addNode(cell, row);

          unsigned y[8];
          y[0] = x[0] + ((lst == 0 && row == 0) ? 1 : 0);
          y[1] = x[1] + ((lst == 0 && row == 1) ? 1 : 0);
          y[2] = x[2] + ((lst == 1 && row == 0) ? 1 : 0);
          y[3] = x[3] + ((lst == 1 && row == 1) ? 1 : 0);
          y[4] = x[4] + ((lst == 2 && row == 0) ? 1 : 0);
          y[5] = x[5] + ((lst == 2 && row == 1) ? 1 : 0);
          y[6] = x[6] + ((lst == 3 && row == 0) ? 1 : 0);
          y[7] = x[7] + ((lst == 3 && row == 1) ? 1 : 0);

          key = y[0];
          for (i = 1; i != 8; i++) {
            key += (y[i] * off[i]);
          }

          if (all_entries[key] == 0) {
            // abkfatal(!in_open[key] && !in_clos[key], "DP2D internal error") ;

            all_entries[key] = getTableEntryFromPool();
            abkassert(all_entries[key] != 0 && all_entries[key]->st != 0,
                      "error: ask Andrew Kennings ");
            visited.push_back(key);

            for (i = 0; i != 8; i++) {
              all_entries[key]->x[i] = y[i];
            }
            all_entries[key]->st->copy(
                init);  // nb: make sure we're clean. need this here???
          } else {
            abkfatal(all_entries[key] != 0, "DP2D internal error");
          }
          // at this point, 'all_entries[key]' holds a pointer to the child
          // which may get updated, or may stay the same depending on bounding.

          tableEntry* child = all_entries[key];

          // is the child already in the open queue?  if no, then insert into
          // the open queue (via heap).  if yes, compare costs.  if new cost
          // is cheaper, we have a less expensive way to the node so we must
          // update and re-heap.

          // we can bound away from crappy solutions.  the amount of violation
          // is how far we have exceeded the valid row width.

          double curr_viol[2];

          curr_viol[0] = max(0.0, (work1.getRiteEdge(0) - rite_edge[0])) +
                         max(0.0, (left_edge[0] - work1.getLeftEdge(0)));

          curr_viol[1] = max(0.0, (work1.getRiteEdge(1) - rite_edge[1])) +
                         max(0.0, (left_edge[1] - work1.getLeftEdge(1)));

          // never return more violation off the end of a row.

          bool bound_away = false;
          unsigned bound_cause = 0;
          if (max(curr_viol[0], curr_viol[1]) >
              max(best_viol[0], best_viol[1])) {
            bound_cause = 1;
            bound_away = true;
          } else {
            if (_usingInitSoln) {
              if (work1.getCost() >= best_cost) {
                bound_cause = 2;
                bound_away = true;
              }
            } else {
              if (max(curr_viol[0], curr_viol[1]) ==
                  max(best_viol[0], best_viol[1])) {
                if (work1.getCost() >= best_cost) {
                  bound_cause = 3;
                  bound_away = true;
                }
              }
            }
          }

          if (!bound_away) {
            if (work1.getLeftEdge(0) < left_edge[0] ||
                work1.getLeftEdge(1) < left_edge[1]) {
              bound_away = true;
            }
          }

          if (bound_away) {
            // cout<<"    bounding away." << endl ;

            continue;
          }

          // see if insersion into the queue is required.

          if (!in_open[key] && !in_clos[key]) {
            all_entries[key]->parent = first;
            all_entries[key]->st->copy(work1);

            // cout<<"------>new insersion: cost:
            // "<<all_entries[key]->st->getCost()<<endl;

            open_queue.push_back(all_entries[key]);
            push_heap(open_queue.begin(), open_queue.end(),
                      two_dp_table_compare());
            in_open[key] = true;
          } else {
            if (in_open[key]) {
              // abkfatal(in_clos[key] == false, "DP2D internal error");

              // cout<<"------>currently in open..."<<endl;

              // nb: if in the open queue, then the cost stored in the keyed
              //     entry will be it's cost in the heap.

              if (work1.getCost() < all_entries[key]->st->getCost()) {
                // find in the heap.

                vector<tableEntry*>::iterator loc;
                for (loc = open_queue.begin(); loc != open_queue.end(); loc++) {
                  if ((*loc) == all_entries[key]) {
                    break;
                  }
                }
                // abkfatal(loc != open_queue.end(), "DP2D internal error");

                (*loc)->parent = first;
                (*loc)->st->copy(work1);
                make_heap(open_queue.begin(), open_queue.end(),
                          two_dp_table_compare());
              } else {
                // cout<<"------>ignore... current cost is better:
                // "<<all_entries[key]->st->getCost()<<endl;
                ;
              }
            } else {
              // abkfatal(in_clos[key] == true, "DP2D internal error");

              // cout<<"------>currently in closed..."<<endl;

              // currently in closed queue; means children have been expanded...
              // re-insersion of this node into open (if cheaper) will be
              // revisit all children too.

              if (work1.getCost() < all_entries[key]->st->getCost()) {
                all_entries[key]->parent = first;
                all_entries[key]->st->copy(work1);

                // cout<<"------>re-insersion: cost:
                // "<<all_entries[key]->st->getCost()<<endl;

                open_queue.push_back(all_entries[key]);
                push_heap(open_queue.begin(), open_queue.end(),
                          two_dp_table_compare());
                in_open[key] = true;
                in_clos[key] = false;
              } else {
                // cout<<"------>ignore... previous cost better:
                // "<<all_entries[key]->st->getCost()<<endl;
                ;
              }
            }
          }
        }
      }
      // dump_open_queue(open_queue);
    }
  }

  // trace back solution here....

  // cout << "nStates visited: " << visited.size() << endl;

  double retval = DBL_MAX;
  if (best_table_entry) {
    vector<tableEntry*> traceback;

    tableEntry* curr = best_table_entry;
    while (curr) {
      traceback.push_back(curr);
      curr = curr->parent;
    }

    curr = traceback.back();

    vector<unsigned> idx(4, 0);
    vector<vector<unsigned> > soln(2);

    twoRowDPStacks work(_problem, _firstNets);

    for (i = traceback.size(); i > 0;) {
      curr = traceback[--i];

      if (i != 0) {
        // determine which coordinate differs.

        unsigned lst = 0;
        unsigned row = 0;

        unsigned coord;
        for (coord = 0; coord != 8; coord++) {
          if (traceback[i]->x[coord] != traceback[i - 1]->x[coord]) {
            lst = coord / 2;
            row = coord % 2;
            break;
          }
        }

        unsigned cell =
            ((lst == 0 || lst == 1) ? posWidthCells[lst - 0][idx[lst]++]
                                    : negWidthCells[lst - 2][idx[lst]++]);

        soln[row].push_back(cell);

        _curSoln[cell].x = work.getLeftEdge(row);
        _curSoln[cell].y = work.getBottomEdge(row);

        work.addNode(cell, row);
      }
    }
    retval = work.getCost();

    // cout<<"confirmation of cost: "<<work.getCost()<<endl;

    // copy new solution over old.
    for (l = 0; l != 2; l++) {
      _list[l].clear();
      for (i = 0; i < soln[l].size(); i++) _list[l].push_back(soln[l][i]);
    }
  }

  _maxMem->update("Dynamic Programming Two-Row SmallPlacer");
  // cleanup.

  for (i = 0; i < visited.size(); i++) {
    releaseTableEntryToPool(all_entries[visited[i]]);
  }

  return retval;
}

void twoRowDPSmallPlacer::buildConnectivity() {
  const SmallPlacementNetlist& netlist = _problem.getNetlist();

  const unsigned numNets = _problem.getNumNets();
  const unsigned numCells = _problem.getNumCells();

  unsigned numPins = numNets * (numCells + 1);

  _adjNets = new TwoRowEdgeAndOffset[numPins];
  _firstNets = new TwoRowEdgeAndOffset*[numCells + 1];

  unsigned c, n, idx = 0;

  // nets (and offsets) adjacent to each node.
  for (c = 0; c < _problem.getNumCells(); c++) {
    _firstNets[c] = _adjNets + idx;
    for (n = 0; n < _problem.getNumNets(); n++) {
      if (!netlist.getPinOffset(c, n).isEmpty()) {
        _adjNets[idx].edge = n;
        _adjNets[idx].pinOffset = netlist.getPinOffset(c, n);
        idx++;
      }
      /**************
                  double pinOffsetX = netlist.getPinOffset(c,n).x ;
                  double pinOffsetY = netlist.getPinOffset(c,n).y ;

                  if(pinOffsetX != DBL_MAX)    //they are connected...
                  {
                      _adjNets[idx].edge = n ;
                      _adjNets[idx].pinOffsetX = pinOffsetX ;
                      _adjNets[idx].pinOffsetY = pinOffsetY ;
                      idx++ ;
                  }
       **************/
    }
  }
  _firstNets[c] = _adjNets + idx;
}

void twoRowDPSmallPlacer::initialSoln() {
  // what should we do for an initial soln?
  ;
}
