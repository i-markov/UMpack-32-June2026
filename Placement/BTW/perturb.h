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

// Created by Igor Markov, Sept 29 1997, VLSICAD ABKGroup UCLA/CS

#ifndef _PERTURB_H_
#define _PERTURB_H_

#include <fstream>
#include <iostream>

#include "abkcommon.h"
#include "placement.h"
#include "rwalk.h"

class Perturbation {
 public:
  // _PreserveRanks randomly shifts cells relative to each other
  // _PreserveCOG   shifts cells independently
  // _PreserveRanks preserves COG too
  enum Type { _PreserveRanks, _PreserveCOG };

 private:
  Type _type;
  double _maxDeltaX;
  double _maxDeltaY;
  bool _haveSeed;
  unsigned _seed;

  RandomWalk *xWalk;
  RandomWalk *yWalk;

  void _initWalks(unsigned size);

  void _perturbPreservingRanks(Placement &);
  void _perturbPreservingCOG(Placement &);

 public:
  Perturbation(Type type, double maxDeltaX, double maxDeltaY);
  Perturbation(Type type, double maxDeltaX, double maxDeltaY, unsigned seedN);
  ~Perturbation();
  void perturb(Placement &);

  unsigned getSeed() const { return xWalk->getSeed(); }
};

#endif
