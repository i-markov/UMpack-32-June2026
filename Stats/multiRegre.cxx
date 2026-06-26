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

// Created 990322 by Stefanus Mantik
#ifdef _MSC_VER
#pragma warning(disable : 4786)
#endif

#include <Stats/multiRegre.h>

#include <cmath>
#include <limits>

#include "statsChecks.h"

using std::cout;
using std::endl;
double MultipleRegression::determinant(const std::vector<std::vector<double> >& m) {
  std::vector<std::vector<double> > e(m);
  double det = 1;
  unsigned i;
  for (i = 0; i < e.size() - 1; i++) {
    if (e[i][i] == 0) {
      unsigned j, jFound = std::numeric_limits<unsigned>::max();
      for (j = i + 1; j < e.size(); j++)
        if (e[j][i] != 0) {
          jFound = j;
          break;
        }
      if (jFound != std::numeric_limits<unsigned>::max())
        for (j = 0; j < e[i].size(); j++) e[i][j] += e[jFound][j];
      else
        det = 0;
    }
    if (det == 0) break;
    for (unsigned j = i + 1; j < e.size(); j++) {
      double xm = e[j][i] / e[i][i];
      for (unsigned k = i; k < e[j].size(); k++)
        e[j][k] = e[j][k] - xm * e[i][k];
    }
  }
  if (det == 0) return det;
  for (i = 0; i < e.size(); i++) det *= e[i][i];
  return det;
}

void MultipleRegression::calculateRegression(
    const std::vector<std::vector<double> >& xIn) {
  statsRequire(!xIn.empty(), "MultipleRegression requires at least one predictor");
  statsRequire(!xIn[0].empty(),
               "MultipleRegression requires non-empty predictor vectors");
  const unsigned numPredictors = xIn.size() - 1;
  const unsigned numSamples = xIn[0].size();
  statsRequire(numSamples > 0, "MultipleRegression requires non-empty sample vectors");
  for (unsigned i = 1; i < xIn.size(); ++i)
    statsRequire(xIn[i].size() == numSamples,
                 "MultipleRegression requires equal-length vectors");

  const unsigned dim = numPredictors + 1;  // intercept + predictors
  std::vector<std::vector<double> > normal(dim, std::vector<double>(dim, 0.0));
  std::vector<double> rhs(dim, 0.0);

  for (unsigned sample = 0; sample < numSamples; ++sample) {
    std::vector<double> row(dim, 1.0);
    for (unsigned predictor = 0; predictor < numPredictors; ++predictor)
      row[predictor + 1] = xIn[predictor][sample];

    const double y = xIn.back()[sample];
    for (unsigned i = 0; i < dim; ++i) {
      rhs[i] += row[i] * y;
      for (unsigned j = 0; j < dim; ++j) normal[i][j] += row[i] * row[j];
    }
  }

  // Solve the normal equations with partial pivoting.
  std::vector<std::vector<double> > aug(dim, std::vector<double>(dim + 1, 0.0));
  for (unsigned i = 0; i < dim; ++i) {
    for (unsigned j = 0; j < dim; ++j) aug[i][j] = normal[i][j];
    aug[i][dim] = rhs[i];
  }

  for (unsigned pivot = 0; pivot < dim; ++pivot) {
    unsigned best = pivot;
    double bestAbs = std::fabs(aug[pivot][pivot]);
    for (unsigned row = pivot + 1; row < dim; ++row) {
      double curAbs = std::fabs(aug[row][pivot]);
      if (curAbs > bestAbs) {
        bestAbs = curAbs;
        best = row;
      }
    }
    statsRequire(bestAbs != 0.0, "Variables are not independent");
    if (best != pivot) {
      for (unsigned col = pivot; col <= dim; ++col) {
        double tmp = aug[pivot][col];
        aug[pivot][col] = aug[best][col];
        aug[best][col] = tmp;
      }
    }

    const double pivotValue = aug[pivot][pivot];
    for (unsigned col = pivot; col <= dim; ++col) aug[pivot][col] /= pivotValue;

    for (unsigned row = 0; row < dim; ++row) {
      if (row == pivot) continue;
      const double factor = aug[row][pivot];
      if (factor == 0.0) continue;
      for (unsigned col = pivot; col <= dim; ++col)
        aug[row][col] -= factor * aug[pivot][col];
    }
  }

  _c = aug[0][dim];
  for (unsigned predictor = 0; predictor < numPredictors; ++predictor)
    _k[predictor] = aug[predictor + 1][dim];

  double sse = 0.0;
  double sst = 0.0;
  double ySum = 0.0;
  for (unsigned sample = 0; sample < numSamples; ++sample) ySum += xIn.back()[sample];
  const double yMean = ySum / numSamples;
  for (unsigned sample = 0; sample < numSamples; ++sample) {
    double predicted = _c;
    for (unsigned predictor = 0; predictor < numPredictors; ++predictor)
      predicted += _k[predictor] * xIn[predictor][sample];
    const double residual = xIn.back()[sample] - predicted;
    sse += residual * residual;
    const double centered = xIn.back()[sample] - yMean;
    sst += centered * centered;
  }
  const double deter = (sst == 0.0) ? 0.0 : 1.0 - sse / sst;
  const double dvest = (numSamples <= dim) ? 0.0 : std::sqrt(sse / (numSamples - dim));
  cout << "determination coeff = " << deter << endl;
  cout << "STDEV of estimative = " << dvest << endl;
}

MultipleRegression::MultipleRegression(const std::vector<std::vector<double> >& x,
                                       const std::vector<double>& y)
    : _k(x.size(), 0.0) {
  statsRequire(!x.empty(), "MultipleRegression requires at least one predictor");
  statsRequire(x[0].size() == y.size(),
               "Attempt to compute linear regression for vectors of different sizes");
  statsRequire(!y.empty(),
               "Attempt to compute linear regression for zero-sized vectors");
  for (unsigned i = 1; i < x.size(); ++i)
    statsRequire(x[i].size() == y.size(),
                 "All predictor vectors must match the response size");

  /*
    unsigned i, size=x.size();
    double   sx = 0, sy = 0;

    for (i=0;i!=size; i++)
    {
      sx += x[i];
      sy += y[i];
    }

    double sxoss = sx/size;
    double st2 = 0;
    for(i=0;i!=size;i++)
    {
      double  t = x[i] - sxoss;
      st2 += t * t;
      _k    += t*y[i];
    }

    _k /= st2;
    _c = (sy-sx*_k)/size;
  */
  std::vector<std::vector<double> > r(x);
  r.push_back(y);
  calculateRegression(r);
}
