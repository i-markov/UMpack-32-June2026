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

// Created: September 22, 1997 by Igor Markov
// CHANGES
#ifdef _MSC_VER
#pragma warning(disable : 4786)
#endif

#include "trivStats.h"

#include <cmath>
#include <iostream>
#include <limits>

#include "statsChecks.h"

using std::endl;
using std::ostream;

Correlation::Correlation(const std::vector<double>& x1,
                         const std::vector<double>& x2) {
  statsRequire(x1.size() == x2.size(),
               "Attempted to compute correlation of vectors of different sizes");
  statsRequire(!x1.empty(),
               "Attempted to compute correlation of zero-sized vectors");
  unsigned k, size = x1.size();
  double x1sum = 0.0, x2sum = 0.0;
  for (k = 0; k != size; k++) {
    x1sum += x1[k];
    x2sum += x2[k];
  }
  double x1ave = x1sum / size, x2ave = x2sum / size;

  double dx1dx1 = 0.0, dx2dx2 = 0.0, dx1dx2 = 0.0;
  for (k = 0; k != size; k++) {
    double dx1 = x1[k] - x1ave, dx2 = x2[k] - x2ave;
    dx1dx1 += dx1 * dx1;
    dx2dx2 += dx2 * dx2;
    dx1dx2 += dx1 * dx2;
  }

  _corr = dx1dx2 / sqrt(dx1dx1 * dx2dx2);
}

Correlation::Correlation(const std::vector<unsigned>& x1,
                         const std::vector<unsigned>& x2) {
  statsRequire(x1.size() == x2.size(),
               "Attempted to compute correlation of vectors of different sizes");
  statsRequire(!x1.empty(),
               "Attempted to compute correlation of zero-sized vectors");
  unsigned k, size = x1.size();
  double x1sum = 0.0, x2sum = 0.0;
  for (k = 0; k != size; k++) {
    x1sum += x1[k];
    x2sum += x2[k];
  }
  double x1ave = x1sum / size, x2ave = x2sum / size;

  double dx1dx1 = 0.0, dx2dx2 = 0.0, dx1dx2 = 0.0;
  for (k = 0; k != size; k++) {
    double dx1 = x1[k] - x1ave, dx2 = x2[k] - x2ave;
    dx1dx1 += dx1 * dx1;
    dx2dx2 += dx2 * dx2;
    dx1dx2 += dx1 * dx2;
  }

  _corr = dx1dx2 / sqrt(dx1dx1 * dx2dx2);
}

TrivialStats::TrivialStats(const std::vector<float>& data)
    : _tot(0.0), _num(data.size()) {
  statsWarn(_num != 0, "Can't create TrivialStats from a zero-sized vector");
  if (_num == 0) {
    _avg = _tot = _max = _min = std::numeric_limits<double>::max();
    return;
  }
  _max = _min = data[0];
  for (unsigned k = 0; k != _num; k++) {
    double cur(data[k]);
    _tot += cur;
    _max = std::max(_max, cur);
    _min = std::min(_min, cur);
  }
  _avg = (_tot + 0.0) /
         _num;  // 0.0 needed in case we reuse this code for int types
}

TrivialStats::TrivialStats(const std::vector<double>& data)
    : _tot(0.0), _num(data.size()) {
  statsWarn(_num != 0, "Can't create TrivialStats from a zero-sized vector");
  if (_num == 0) {
    _avg = _tot = _max = _min = std::numeric_limits<double>::max();
    return;
  }
  _max = _min = data[0];
  for (unsigned k = 0; k != _num; k++) {
    double cur(data[k]);
    _tot += cur;
    _max = std::max(_max, cur);
    _min = std::min(_min, cur);
  }
  _avg = (_tot + 0.0) /
         _num;  // 0.0 needed in case we reuse this code for int types
}

TrivialStats::TrivialStats(const std::vector<unsigned>& data)
    : _tot(0.0), _num(data.size()) {
  statsRequire(_num != 0,
               "Can't create TrivialStats from a zero-sized vector");
  _max = _min = data[0];
  for (unsigned k = 0; k != _num; k++) {
    double cur(data[k]);
    _tot += cur;
    _max = std::max(_max, cur);
    _min = std::min(_min, cur);
  }
  _avg = (_tot + 0.0) /
         _num;  // 0.0 needed in case we reuse this code for int types
}

ostream& operator<<(ostream& out, const TrivialStats& stats) {
  out << " Max: " << stats.getMax() << "  Min: " << stats.getMin()
      << "  Avg: " << stats.getAvg() << endl;
  return out;
}

TrivialStatsWithStdDev::TrivialStatsWithStdDev(const std::vector<float>& data)
    : TrivialStats(data) {
  float totDev = 0.0;
  for (unsigned k = 0; k != _num; k++) totDev += (_avg - data[k]) * (_avg - data[k]);
  _stdDev = sqrt(totDev / _num);
}

TrivialStatsWithStdDev::TrivialStatsWithStdDev(const std::vector<double>& data)
    : TrivialStats(data) {
  double totDev = 0.0;
  for (unsigned k = 0; k != _num; k++) totDev += (_avg - data[k]) * (_avg - data[k]);
  _stdDev = sqrt(totDev / _num);
}

TrivialStatsWithStdDev::TrivialStatsWithStdDev(const std::vector<unsigned>& data)
    : TrivialStats(data) {
  double totDev = 0.0;
  for (unsigned k = 0; k != _num; k++) totDev += (_avg - data[k]) * (_avg - data[k]);
  _stdDev = sqrt(totDev / _num);
}

ostream& operator<<(ostream& out, const TrivialStatsWithStdDev& stats) {
  out << " Max: " << stats.getMax() << "  Min: " << stats.getMin()
      << "  Avg: " << stats.getAvg() << "  StdDev: " << stats.getStdDev()
      << endl;
  return out;
}
