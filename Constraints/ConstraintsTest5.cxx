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

// Fixed-type constraint application driver.
// It prepares a placement, copies the fixed-type constraints onto a larger
// target placement, enforces them, and prints the resulting constraint state
// and fixed bounding box.

#include <iostream>

#include "Placement/subPlacement.h"
#include "allConstr.h"
#include "constraints.h"

using std::cout;
using std::endl;

int main() {
  Constraints ctrs;

  Placement pl(5);
  pl[0] = pl[1] = pl[2] = pl[3] = pl[4] = Point(0.4444, 0.4444);
  pl[2] = Point(0.333, 0.333);
  pl[0] = Point(0.999, 0.999);
  //  cout << "------ Original placement " << endl;
  //  cout << pl;

  Mapping mp(3, 5);
  mp[0] = 2, mp[1] = 0, mp[2] = 4;
  Mapping mp3(2, 5);
  mp3[0] = 2, mp3[1] = 4;
  //  cout << "------ Original subset " << endl << mp << endl;

  {
    // NOTE: simplified ctors have been added to FixedTypeConstraints lately.
    //       so that you do not have to construct a SubPlacement
    //       If you understand the examples here, look up the simplified ctors
    //       in fxTypeConstr.h  -ILM

    SubPlacement subPl(pl, mp3);

    Placement fixedLocs(2);
    fixedLocs[0] = Point(0, 0);
    fixedLocs[1] = Point(1, 1);

    FixedGroupConstraint fgc(subPl, fixedLocs);
    ctrs.add(fgc);
  }
  {
    SubPlacement subPl(pl, mp);

    Placement fixedLocs(3);
    fixedLocs[0] = Point(0, 0);
    fixedLocs[1] = Point(1, 1);
    fixedLocs[2] = Point(2, 2);

    TetheredGroupConstraint fgc(subPl, fixedLocs, 0.1);
    ctrs.add(fgc);
  }
  {
    SubPlacement subPl(pl, mp);
    BBox bb(0.25, 0.25, 0.35, 0.35);
    RectRegionConstraint rrc(subPl, bb);
    ctrs.add(rrc);
  }
  {
    SubPlacement subPl(pl, mp);
    EqualXConstraint cgc(subPl);
    ctrs.add(cgc);
  }
  {
    SubPlacement subPl(pl, mp);
    EqualYConstraint cgc(subPl);
    ctrs.add(cgc);
  }
  {
    SubPlacement subPl(pl, mp);
    SoftGroupConstraint cgc(subPl, 1, 1, 1);
    ctrs.add(cgc);
  }

  // Additional FixedGroup-type constraints
  Mapping mp1(2, 5);
  mp1[0] = 1, mp1[1] = 4;
  Mapping mp2(1, 5);
  mp2[0] = 3;
  {
    SubPlacement subPl(pl, mp1);
    Placement fixedLocs(2);
    fixedLocs[0] = Point(0.1, 0.1);
    fixedLocs[1] = Point(1, 1);

    FixedGroupConstraint fgc(subPl, fixedLocs);
    ctrs.add(fgc);
  }
  {
    SubPlacement subPl(pl, mp2);
    Placement fixedLocs(1);
    fixedLocs[0] = Point(0.1, 0.1);

    FixedXGroupConstraint fgc(subPl, fixedLocs);
    ctrs.add(fgc);
  }

  ctrs.compactify();
  cout << " Compactified: " << ctrs << endl;

  if (ctrs.hasFGC()) {
    cout << " ===== The FixedGroup constraint ====== " << endl;
    cout << ctrs.getFGC().getFixedLocs() << endl;
    cout << ctrs.getFGC().getMapping() << endl;
  }

  Placement pl1(10);
  Constraints ctrs1(ctrs, pl1, Constraints::FIXEDTYPE);

  ctrs1.compactify();
  cout << " New constraints: " << ctrs1 << endl;
  cout << " Placement size : " << ctrs1.getPlacement().getSize() << endl;

  if (ctrs1.hasFGC()) {
    cout << " ===== The FixedGroup constraint ====== " << endl;
    cout << ctrs1.getFGC().getFixedLocs() << endl;
    cout << ctrs1.getFGC().getMapping() << endl;
  }

  cout << " Fixed BBox : " << ctrs1.getFixedBBox() << endl;

  ctrs1.enforceAll();

  Constraints ctrs2(ctrs, pl1, Constraints::ALL);
  cout << ctrs2;
}
