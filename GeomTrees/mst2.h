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

/*                                                                    */
/*            (c) Copyright, 1997 by Professor Gabriel Robins         */
/*                                                                    */
/*      Department of Computer Science, University of Virginia        */
/*          Charlottesville, VA 22903-2442  (804) 982-2207            */
/*   robins@cs.virginia.edu    http://www.cs.virginia.edu/~robins/    */
/*                                                                    */
/*   This code may be freely used for all non-commercial purposes.    */
/*   All copies/portions of this code must contain this header.       */
/*                                                                    */

// This code is a C++ implementation of MST2 from I1Steiner package
// The code has been modified by Stefanus Mantik from C to C++

#include "Geoms/plGeom.h"

class MST2 {
 public:
  int P1, P2;
  double weight;

  I1StEdge() : P1(0), P2(0), weight(DBL_MAX) {}
};

class I1StNode {
 public:
  int parent, id;
  double weight;

  I1StNode() : parent(-1), id(-1), weight(DBL_MAX) {}
};

/**********************************************************************/

class I1Steiner {
 public:
  unsigned nPoints;

 private:
  vector<Point> pointset, SP_candidates;
  I1StNode *Tree, *old_Tree;
  vector<I1StEdge> Neighbour;
  vector<unsigned> degree_count;

  vector<double> Steiner_savings;
  unsigned id_num;
  bool _steinerComputed;

  void init_vars();
  void qsort_candidates(int, int);
  void closest_neighbour();
  void find_neighbours(int);
  void update_Tree(int, int);
  void move_up(int &, int &);
  int reset_po(int, int);
  double linear_MST(int);
  void remove_invalid_SPs(unsigned, unsigned &);

  void generate_MST3();
  double generate_steiner5();

 public:
  double mst_cost;

  I1Steiner(unsigned n_points);
  I1Steiner(const vector<Point> &points);
  ~I1Steiner();

  double computeSteiner() {
    if (_steinerComputed) return mst_cost;
    generate_MST3();
    return generate_steiner5();
  }

  void print_points();
  void add_point(unsigned index, double X, double Y);
};

/************************************************************************/

#ifndef _SWAP_
#define _SWAP_
template <class T>
inline void SWAP(T &a, T &b) {
  T temp;
  temp = a;
  a = b;
  b = temp;
}
#endif
