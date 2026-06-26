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

// Created by Igor Markov, July 1997, VLSICAD ABKGroup UCLA/CS

#include "perturb.h"

#include <fstream>
#include <iostream>

#include "abkcommon.h"
#include "placement.h"
#include "rwalk.h"

Perturbation::Perturbation(Type type, double maxDeltaX, double maxDeltaY)
    : _type(type),
      _maxDeltaX(maxDeltaX),
      _maxDeltaY(maxDeltaY),
      _haveSeed(false),
      xWalk(NULL),
      yWalk(NULL) {}

Perturbation::Perturbation(Type type, double maxDeltaX, double maxDeltaY,
                           unsigned seedN)
    : _type(type),
      _maxDeltaX(maxDeltaX),
      _maxDeltaY(maxDeltaY),
      _haveSeed(true),
      _seed(seedN),
      xWalk(NULL),
      yWalk(NULL) {}

void Perturbation::_initWalks(unsigned size) {
  if (xWalk) delete xWalk;
  if (yWalk) delete yWalk;
  if (_haveSeed) {
    xWalk = new RandomWalk(size, _maxDeltaX, 0.0, _seed);
    yWalk = new RandomWalk(size, _maxDeltaY, 0.0,
                           xWalk->getSeed() * xWalk->getSeed());
  } else {
    xWalk = new RandomWalk(size, _maxDeltaX, 0.0);
    yWalk = new RandomWalk(size, _maxDeltaY, 0.0,
                           xWalk->getSeed() * xWalk->getSeed());
  }
}

Perturbation::~Perturbation() {
  if (xWalk) delete xWalk;
  if (yWalk) delete yWalk;
}

void Perturbation::perturb(Placement& pl) {
  abkfatal(pl.getSize() != 0, " Attempting to perturb a Placment of size 0 ");
  _initWalks(pl.getSize());
  switch (_type) {
    case _PreserveRanks:
      _perturbPreservingRanks(pl);
      break;
    case _PreserveCOG:
      _perturbPreservingCOG(pl);
      break;
    default:
      abkfatal(0, " Unknown perturbation type ");
  }
}

void Perturbation::_perturbPreservingRanks(Placement& pl) {
  unsigned k, size = pl.getSize();
  vector<double> vecX(size), vecY(size);
  for (k = 0; k != size; k++) {
    vecX[k] = pl[k].x;
    vecY[k] = pl[k].y;
  }

  Permutation xIdx;
  Permutation yIdx;

  Permutation xInv(vecX);
  Permutation yInv(vecY);

  xInv.getInverse(xIdx);
  yInv.getInverse(yIdx);

  double epsilon = 1e-4;
  double xLag = 0.0, yLag = 0.0;
  for (k = 1; k != size; k++) {
    double newX = pl[xIdx[k]].x + xWalk->getLoc(k) - xLag;
    double newY = pl[yIdx[k]].y + yWalk->getLoc(k) - yLag;
    if (newX < pl[xIdx[k - 1]].x) {
      xLag = pl[xIdx[k - 1]].x + epsilon - newX;
      newX = pl[xIdx[k - 1]].x + epsilon;
    } else
      xLag = 0.0;
    if (newY < pl[yIdx[k - 1]].y) {
      yLag = pl[yIdx[k - 1]].y + epsilon - newY;
      newY = pl[yIdx[k - 1]].y + epsilon;
    } else
      yLag = 0.0;

    pl[xIdx[k]].x = newX;
    pl[yIdx[k]].y = newY;
  }
}

void Perturbation::_perturbPreservingCOG(Placement& pl) {
  unsigned size = pl.getSize();
  for (unsigned k = 0; k != size; k++) {
    pl[k].x += xWalk->getStep(k);
    pl[k].y += yWalk->getStep(k);
  }
}
