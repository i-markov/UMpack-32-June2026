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

// Created: July 5, 1997 by Igor Markov  VLSICAD ABKGroup UCLA/CS
#ifdef _MSC_VER
#pragma warning(disable : 4786)
#endif

#include "compar.h"

#include "placementUtils.h"

namespace {

struct PlacementPointAccessor {
  const Placement& pl;
  PlacementPointAccessor(const Placement& p) : pl(p) {}
  const Point& operator()(unsigned i) const { return pl[i]; }
};

template <class Accumulator>
double samplePairMetric(Placement& pl1, Placement& pl2, Accumulator acc) {
  RandomUnsigned randIdx(0, pl1.getSize());
  double value = acc.initial(pl1.getSize());
  for (unsigned i = 0; i < 5 * pl1.getSize(); i++) {
    unsigned idxA = randIdx;
    unsigned idxB = randIdx;
    acc.step(value, mDist(pl1[idxA], pl1[idxB]), mDist(pl2[idxA], pl2[idxB]));
  }
  return acc.finish(value, pl1.getSize());
}

struct AvgElongationAccumulator {
  double initial(unsigned) const { return 0.0; }
  void step(double& value, double dist1, double dist2) const {
    value += fabs(dist2 - dist1);
  }
  double finish(double value, unsigned nPts) const {
    return value / (5 * nPts);
  }
};

struct MaxElongationAccumulator {
  double initial(unsigned) const { return 0.0; }
  void step(double& value, double dist1, double dist2) const {
    const double delta = fabs(dist2 - dist1);
    if (value < delta) value = delta;
  }
  double finish(double value, unsigned) const { return value; }
};

struct AvgStretchAccumulator {
  double initial(unsigned) const { return 0.0; }
  void step(double& value, double dist1, double dist2) const {
    if (dist1 > 1e-15 && dist2 > 1e-15) {
      double rat = dist2 / dist1;
      value += (rat > 1) ? rat : 1 / rat;
    }
  }
  double finish(double value, unsigned nPts) const {
    return value / (5 * nPts);
  }
};

struct MaxStretchAccumulator {
  double initial(unsigned) const { return 0.0; }
  void step(double& value, double dist1, double dist2) const {
    if (dist1 > 1e-15 && dist2 > 1e-15) {
      double rat = dist2 / dist1;
      if (rat > 1) rat = 1 / rat;
      if (value < rat) value = rat;
    }
  }
  double finish(double value, unsigned) const { return value; }
};

}  // namespace

double MaxMDiff ::operator()(Placement& pl1, Placement& pl2) {
  abkfatal(pl1.getSize() == pl2.getSize(),
           "Can\'t compare placements of different sizes");
  return placement_utils::linfDistance(
      pl1.getSize(), PlacementPointAccessor(pl1), PlacementPointAccessor(pl2));
}

double AvgMDiff::operator()(Placement& pl1, Placement& pl2) {
  abkfatal(pl1.getSize() == pl2.getSize(),
           "Can\'t compare placements of different sizes");
  return placement_utils::averageDistance(
      pl1.getSize(), PlacementPointAccessor(pl1), PlacementPointAccessor(pl2));
}

double AvgElongation::operator()(Placement& pl1, Placement& pl2) {
  abkfatal(pl1.getSize() == pl2.getSize(),
           "Can\'t compare placements of different sizes");
  return samplePairMetric(pl1, pl2, AvgElongationAccumulator());
}

double MaxElongation::operator()(Placement& pl1, Placement& pl2) {
  abkfatal(pl1.getSize() == pl2.getSize(),
           "Can\'t compare placements of different sizes");
  return samplePairMetric(pl1, pl2, MaxElongationAccumulator());
}

double AvgStretch::operator()(Placement& pl1, Placement& pl2) {
  abkfatal(pl1.getSize() == pl2.getSize(),
           "Can\'t compare placements of different sizes");
  return samplePairMetric(pl1, pl2, AvgStretchAccumulator());
}

double MaxStretch::operator()(Placement& pl1, Placement& pl2) {
  abkfatal(pl1.getSize() == pl2.getSize(),
           "Can\'t compare placements of different sizes");
  return samplePairMetric(pl1, pl2, MaxStretchAccumulator());
}
