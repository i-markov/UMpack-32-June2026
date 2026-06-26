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
***  EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
***  MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
***  IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
***  CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
***  TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
***  SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
***
***************************************************************************/

// Created:    Igor Markov,  VLSI CAD ABKGROUP UCLA  Sept 29 1997

// 971113 mro unified seeded, unseeded ctors in RandomWalk

#ifndef _RWALK_H_
#define _RWALK_H_

#include "abkcommon.h"

class RandomWalk;
typedef RandomWalk RWalk;

class RandomWalk
// One-dimensional fixed length random walk starting at 0.0
// and returning to a predefined location (default 0.0)
// with steps limited by _maxStep in each direction
//
// if needed, I can increase the quality of "randomness" by an additional
// random shuffle of random steps
{
  double _maxStep;
  double _endsAt;
  vector<double> _locations;
  vector<double> _steps;
  RandomDouble _randouble;

  void _generate();

 public:
  RandomWalk(unsigned size, double maxStep, double endsAt = 0.0,
             unsigned seedN = UINT_MAX);

  // to "reconstruct events"
  unsigned getSeed() const { return _randouble.getSeed(); }

  // synonyms
  unsigned getSize() const { return _steps.size(); }
  unsigned getLength() const { return _steps.size(); }

  double getStep(unsigned idx) const { return _steps[idx]; }
  double getLoc(unsigned idx) const { return _locations[idx]; }

  const vector<double>& getSteps() const { return _steps; }
  const vector<double>& getLocs() const { return _locations; }
};

ostream& operator<<(ostream&, const RandomWalk&);

#endif
