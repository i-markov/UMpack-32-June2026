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

#include "bfsGen.h"

#include <algorithm>
#include <deque>
#include <vector>

#include "HGraph/hgFixed.h"
#include "Stats/stats.h"
#include "ordering.h"

using std::deque;
using std::shuffle;
using std::vector;

// SBGWeightBfs: implemented by DAP
// Produces a Breadth First Search based ordering on the vertices of the
// hypergraph, and then partitions according to tolerances.

vector<double> SBGWeightBfs::generateSoln(Partitioning& curPart) {
  unsigned numPartitions = _problem.getNumPartitions();
  const HGraphFixed& hgraph = _problem.getHGraph();
  abkfatal(numPartitions == 2,
           "this solution generator (Bfs) works only for"
           " 2-way partitioning problems");

  abkfatal(hgraph.getNumNodes() >= numPartitions,
           "not enough nodes / partition");

  vector<double> areaToMin(numPartitions), areaToFill(numPartitions),
      targetAreas(numPartitions), availableArea(numPartitions),
      chancesOfGoingToPart(numPartitions);

  const std::vector<std::vector<double> >& maxCaps =
      _problem.getMaxCapacities();
  const std::vector<std::vector<double> >& minCaps =
      _problem.getMinCapacities();
  const std::vector<std::vector<double> >& caps = _problem.getCapacities();

  unsigned k;
  for (k = 0; k < numPartitions; k++) {
    availableArea[k] = maxCaps[k][0];
    areaToMin[k] = minCaps[k][0];
    areaToFill[k] = caps[k][0];
    targetAreas[k] = caps[k][0];
  }

  int v_ct = hgraph.getNumNodes(), e_ct = hgraph.getNumEdges();
  std::vector<short> vvisited(v_ct, false);
  std::vector<short> evisited(e_ct, false);
  int current;
  int last = 0;
  bool done = false, front_turn = true;
  ordering ord(v_ct);
  deque<int> frontQ, backQ;
  if (hgraph.getNumTerminals() == 2) {
    frontQ.push_back((*hgraph.terminalsBegin())->getIndex());
    vvisited[((*hgraph.terminalsBegin())->getIndex())] = true;
    backQ.push_back((*(hgraph.terminalsEnd() - 1))->getIndex());
    vvisited[((*(hgraph.terminalsEnd() - 1))->getIndex())] = true;

  } else if (hgraph.getNumTerminals() == 1) {
    frontQ.push_back((*hgraph.terminalsBegin())->getIndex());
    vvisited[((*hgraph.terminalsBegin())->getIndex())] = true;
  } else {
    RandomUnsigned r(0, v_ct);
    unsigned f1 = r, b1 = r;
    frontQ.push_back(f1);
    backQ.push_back(b1);
    vvisited[f1] = true;
    vvisited[b1] = true;
  }
  for (int i = (hgraph.terminalsEnd() - hgraph.terminalsBegin()); i < v_ct; i++)
    // treat all nodes as terminals which are only in one part
    if (!curPart[i].isInPart(0) && curPart[i].isInPart(1)) {
      frontQ.push_back(i);
      vvisited[i] = true;
    } else if (curPart[i].isInPart(0) && !curPart[i].isInPart(1)) {
      backQ.push_back(i);
      vvisited[i] = true;
    } else if (!curPart[i].isInPart(0) && !curPart[i].isInPart(1)) {
      abkfatal(0, "No possible solution!");
    }

  int fronti = 0, backi = v_ct - 1;
  // initialize queues
  RandomUnsigned r(0, UINT_MAX);
  do {
    while (frontQ.size() > 0 || backQ.size() > 0) {
      if (front_turn && frontQ.size() == 0) front_turn = false;
      if (!front_turn && backQ.size() == 0) front_turn = true;
      if (front_turn) {
        current = frontQ.front();
        frontQ.pop_front();
      } else {
        current = backQ.front();
        backQ.pop_front();
      }
      if (!vvisited[current])
        vvisited[current] = true;  // get this out of here if possible
      int front_wave_i = -1;
      if (front_turn) {
        ord.set(fronti++, current);
        front_wave_i = frontQ.size();
      } else {
        ord.set(backi--, current);
        front_wave_i = backQ.size();
      }

      const HGFNode& curr_n = hgraph.getNodeByIdx(current);
      ctainerHGFEdgesLocal E(curr_n.edgesBegin(), curr_n.edgesEnd());
      sort(E.begin(), E.end());
      for (itHGFEdgeLocal it1 = E.begin(); it1 != E.end(); ++it1) {
        if (!evisited[(*it1)->getIndex()] &&
            (*it1)->getDegree() <= v_ct * 0.10) {
          evisited[(*it1)->getIndex()] = true;
          for (itHGFNodeLocal it2 = (*it1)->nodesBegin();
               it2 != (*it1)->nodesEnd(); ++it2)
            if (!vvisited[(*it2)->getIndex()]) {
              vvisited[(*it2)->getIndex()] = true;
              if (front_turn)
                frontQ.push_back((*it2)->getIndex());
              else
                backQ.push_back((*it2)->getIndex());
            }
        }
      }

      deque<int>::iterator wave_it;
      if (front_turn) {
        wave_it = frontQ.begin() + front_wave_i;
        shuffle(wave_it, frontQ.end(), makeShuffleRng(r));
      } else {
        wave_it = backQ.begin() + front_wave_i;
        shuffle(wave_it, backQ.end(), makeShuffleRng(r));
      }
      front_turn = (r & 1u);  // choose a random next Q
    }

    int j;
    for (j = last; j < v_ct; ++j)
      if (!vvisited[j]) {
        frontQ.push_back(j);
        vvisited[j] = true;
        break;
      }
    last = j + 1;
    if (last >= v_ct) done = true;

  } while (!done);

  vector<int> cuts(v_ct, 0);
  int minp = INT_MAX, maxp = -1;
  int minval = INT_MAX, maxval = -1;
  for (int e = 0; e < e_ct; ++e) {
    const HGFEdge& ed = hgraph.getEdgeByIdx(e);
    for (itHGFNodeLocal it = ed.nodesBegin(); it != ed.nodesEnd(); ++it) {
      if (ord.idx((*it)->getIndex()) < minp) {
        minval = (*it)->getIndex();
        minp = ord.idx((*it)->getIndex());
      }
      if (ord.idx((*it)->getIndex()) > maxp) {
        maxval = (*it)->getIndex();
        maxp = ord.idx((*it)->getIndex());
      }
    }
    cuts[minp]++;
    cuts[maxp]--;
  }
  // end paste

  double start = areaToMin[0];
  double end = availableArea[0];
  double mid = (start + end) / 2;
  double dist_from_mid_tiebreaker = mid - start + 1.0;
  double range = 0.0;

  int icut = 0, mincut = INT_MAX;
  int it = 0, minloc = -1;
  while (range < start) {
    icut += cuts[it];
    const HGFNode& node = hgraph.getNodeByIdx(ord.loc(it++));
    double nodeWt = hgraph.getWeight(node.getIndex());
    range += nodeWt;
  }
  do {
    double currdist = fabs(range - mid);
    if ((icut < mincut) ||
        (icut == mincut && currdist < dist_from_mid_tiebreaker)) {
      mincut = icut;
      minloc = it;
      dist_from_mid_tiebreaker = currdist;
    }
    icut += cuts[it];
    const HGFNode& node = hgraph.getNodeByIdx(ord.loc(it++));
    double nodeWt = hgraph.getWeight(node.getIndex());
    range += nodeWt;
  } while (range < end);

  vector<double> partSizes(numPartitions);

  int o = ord.size();
  for (int i = 0; i < o; i++)
    if (ord.idx(i) < minloc) {
      curPart[i].setToPart(0);
      const HGFNode& node = hgraph.getNodeByIdx(ord.idx(i));
      double nodeWt = hgraph.getWeight(node.getIndex());
      partSizes[0] += nodeWt;
    } else {
      curPart[i].setToPart(1);
      const HGFNode& node = hgraph.getNodeByIdx(ord.idx(i));
      double nodeWt = hgraph.getWeight(node.getIndex());
      partSizes[1] += nodeWt;
    }

  return partSizes;
}
