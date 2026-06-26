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

// 971113 mro removed unseeded ctor for RandomWalk (use default value)

#include "rwalk.h"

#include "abkcommon.h"

RandomWalk::RandomWalk(unsigned size, double maxStep, double endsAt,
                       unsigned seedN)
    : _locations(size),
      _steps(size),
      _maxStep(maxStep),
      _endsAt(endsAt),
      _randouble(0, 1, seedN) {
  abkfatal(size != 0, "Can't create a RandomWalk with 0 steps");
  abkfatal(abs(_endsAt) / size <= _maxStep,
           " Impossible random walk : step too small ");
  _generate();
}
void RandomWalk::_generate() {
  double sum = 0.0;
  unsigned k = 0, size = _steps.size();
  while (k < size - 1)  // note: k can decrease
  {
    double maxReserve = _maxStep * (size - k - 1) - (sum - _endsAt);
    double minReserve = -_maxStep * (size - k - 1) - (sum - _endsAt);
    if (maxReserve < -_maxStep || minReserve > _maxStep) {
      abkfatal(k != 0,
               " Random Walk failed to initialize, email to "
               "igor.markov1@gmail.com ");
      k--;
      continue;
    }
    double hi = min(maxReserve, _maxStep);
    double lo = max(minReserve, -_maxStep);

    _locations[k++] = sum += _steps[k] = lo + _randouble * (hi - lo);
  }

  if (sum - _endsAt < -_maxStep || sum - _endsAt > _maxStep)
    abkfatal(
        0,
        " Random Walk failed to initialize, email to igor.markov1@gmail.com ");

  _locations[k] = sum += _steps[k] = -(sum - _endsAt);
}

ostream& operator<<(ostream& out, const RandomWalk& rwk) {
  unsigned size = rwk.getSize();
  out << " Random walk of length " << size << endl;
  for (unsigned k = 0; k != size; k++)
    out << setw(10) << rwk.getStep(k) << setw(10) << rwk.getLoc(k) << "\n";
  return out;
}
