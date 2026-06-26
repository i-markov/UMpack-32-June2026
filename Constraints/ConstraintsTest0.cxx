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

// Constraint regression driver.
// It instantiates each major constraint type on a small placement, prints
// whether the constraint is satisfied, and shows the effect of enforce() for
// the fixed, tethered, and region-based examples.

#include <iostream>

#include "Combi/mapping.h"
#include "Placement/placement.h"
#include "Placement/subPlacement.h"
#include "allConstr.h"
#include "constraints.h"

using std::cout;
using std::endl;
using std::vector;

int main() {
  Placement pl(5);
  pl[0] = pl[1] = pl[2] = pl[3] = pl[4] = Point(0.4444, 0.4444);
  pl[2] = Point(0.333, 0.333);
  pl[0] = Point(0.999, 0.999);
  cout << "------ Original placement " << endl;
  cout << pl;

  Mapping mp(3, 5);
  mp[0] = 2, mp[1] = 0, mp[2] = 4;
  cout << "------ Original subset " << endl << mp << endl;

  {
    // NOTE: simplified ctors have been added to FixedTypeConstraints lately.
    //       If you understand the examples here, look up the simplified syntax
    //       in fxTypeConstr.h  -ILM

    cout << "----- Testing FixedGroupConstraint " << endl;
    Placement pl1(pl);
    SubPlacement subPl(pl1, mp);

    Placement fixedLocs(3);
    fixedLocs[0] = Point(0, 0);
    fixedLocs[1] = Point(1, 1);
    fixedLocs[2] = Point(2, 2);

    FixedGroupConstraint fgc(subPl, fixedLocs);

    cout << "Is satisfied : " << fgc.isSatisfied() << endl;
    cout << "Penalty      : " << fgc.getPenalty() << endl;
    fgc.enforce();
    cout << "Enforce now \n";
    cout << pl1;
    cout << "Is satisfied : " << fgc.isSatisfied() << endl;
    cout << "Penalty      : " << fgc.getPenalty() << endl;
  }
  {
    cout << "----- Testing FixedXGroupConstraint " << endl;
    Placement pl1(pl);
    SubPlacement subPl(pl1, mp);

    Placement fixedLocs(3);
    fixedLocs[0] = Point(0, 0);
    fixedLocs[1] = Point(1, 1);
    fixedLocs[2] = Point(2, 2);

    FixedXGroupConstraint fgc(subPl, fixedLocs);

    cout << "Is satisfied : " << fgc.isSatisfied() << endl;
    cout << "Penalty      : " << fgc.getPenalty() << endl;
    fgc.enforce();
    cout << "Enforce now \n";
    cout << pl1;
    cout << "Is satisfied : " << fgc.isSatisfied() << endl;
    cout << "Penalty      : " << fgc.getPenalty() << endl;
  }
  {
    cout << "----- Testing FixedYGroupConstraint " << endl;
    Placement pl1(pl);
    SubPlacement subPl(pl1, mp);

    Placement fixedLocs(3);
    fixedLocs[0] = Point(0, 0);
    fixedLocs[1] = Point(1, 1);
    fixedLocs[2] = Point(2, 2);

    FixedYGroupConstraint fgc(subPl, fixedLocs);

    cout << "Is satisfied : " << fgc.isSatisfied() << endl;
    cout << "Penalty      : " << fgc.getPenalty() << endl;
    fgc.enforce();
    cout << "Enforce now \n";
    cout << pl1;
    cout << "Is satisfied : " << fgc.isSatisfied() << endl;
    cout << "Penalty      : " << fgc.getPenalty() << endl;
  }
  {
    cout << "----- Testing TetheredGroupConstraint " << endl;
    Placement pl1(pl);
    SubPlacement subPl(pl1, mp);

    Placement fixedLocs(3);
    fixedLocs[0] = Point(0, 0);
    fixedLocs[1] = Point(1, 1);
    fixedLocs[2] = Point(2, 2);

    TetheredGroupConstraint fgc(subPl, fixedLocs, 0.1);

    cout << "Is satisfied : " << fgc.isSatisfied() << endl;
    cout << "Penalty      : " << fgc.getPenalty() << endl;
    fgc.enforce();
    cout << "Enforce now \n";
    cout << pl1;
    cout << "Is satisfied : " << fgc.isSatisfied() << endl;
    cout << "Penalty      : " << fgc.getPenalty() << endl;
  }
  {
    cout << "----- Testing RectRegionConstraint " << endl;
    Placement pl1(pl);
    SubPlacement subPl(pl1, mp);
    Rectangle bb(0.25, 0.25, 0.35, 0.35);
    RectRegionConstraint rrc(subPl, bb);

    cout << "Is satisfied : " << rrc.isSatisfied() << endl;
    cout << "Penalty      : " << rrc.getPenalty() << endl;
    rrc.enforce();
    cout << "Enforce now \n";
    cout << pl1;
    cout << "Is satisfied : " << rrc.isSatisfied() << endl;
    cout << "Penalty      : " << rrc.getPenalty() << endl;
  }
  {
    cout << "----- Testing COGConstraint " << endl;
    Placement pl1(pl);
    SubPlacement subPl(pl1, mp);
    COGConstraint cgc(subPl, Point(0.5, 0.5));

    cout << "Is satisfied : " << cgc.isSatisfied() << endl;
    cout << "Penalty      : " << cgc.getPenalty() << endl;
    cout << " Readjusting placement\n";
    subPl[0] = Point(0.7, 0.7);
    subPl[1] = Point(0.3, 0.5);
    subPl[2] = Point(0.5, 0.3);
    cout << pl1;
    cout << "Is satisfied : " << cgc.isSatisfied() << endl;
    cout << "Penalty      : " << cgc.getPenalty() << endl;
  }
  {
    cout << "----- Testing EqualXConstraint " << endl;
    Placement pl1(pl);
    SubPlacement subPl(pl1, mp);
    EqualXConstraint cgc(subPl);

    cout << "Is satisfied : " << cgc.isSatisfied() << endl;
    cout << "Penalty      : " << cgc.getPenalty() << endl;
    cout << " Enforce now\n";
    cgc.enforce();
    cout << pl1;
    cout << "Is satisfied : " << cgc.isSatisfied() << endl;
    cout << "Penalty      : " << cgc.getPenalty() << endl;
  }
  {
    cout << "----- Testing EqualYConstraint " << endl;
    Placement pl1(pl);
    SubPlacement subPl(pl1, mp);
    EqualYConstraint cgc(subPl);

    cout << "Is satisfied : " << cgc.isSatisfied() << endl;
    cout << "Penalty      : " << cgc.getPenalty() << endl;
    cout << " Enforce now\n";
    cgc.enforce();
    cout << pl1;
    cout << "Is satisfied : " << cgc.isSatisfied() << endl;
    cout << "Penalty      : " << cgc.getPenalty() << endl;
  }
  {
    cout << "----- Testing StayTogetherConstraint " << endl;
    Placement pl1(pl);
    SubPlacement subPl(pl1, mp);
    Rectangle bbx(0, 0, 0.1, 0.1);
    StayTogetherConstraint cgc(subPl, bbx);

    cout << "Is satisfied : " << cgc.isSatisfied() << endl;
    cout << "Penalty      : " << cgc.getPenalty() << endl;
    cout << " Enforce now\n";
    cgc.enforce();
    cout << pl1;
    cout << "Is satisfied : " << cgc.isSatisfied() << endl;
    cout << "Penalty      : " << cgc.getPenalty() << endl;
  }
  {
    cout << "----- Testing SoftGroupConstraint " << endl;
    Placement pl1(pl);
    SubPlacement subPl(pl1, mp);
    SoftGroupConstraint cgc(subPl, 1, 1, 1);

    pl1[0] = Point(0.5, 0.6);
    pl1[2] = Point(0, 0);
    cout << "Is satisfied : " << cgc.isSatisfied() << endl;
    cout << "Penalty      : " << cgc.getPenalty() << endl;
  }
  {
    cout << "----- Testing MultRectRegionConstraint " << endl;
    Placement pl1(pl);
    pl1[0] = Point(0, 0);
    pl1[2] = Point(0, 0);
    pl1[4] = Point(0.5, 0.6);
    SubPlacement subPl(pl1, mp);
    vector<Rectangle> regions(1);
    regions[0] = BBox(0, 0, 0.5, 0.5);
    MultRectRegionConstraint cgc(subPl, regions);

    cout << "Is satisfied : " << cgc.isSatisfied() << endl;
    cout << "Penalty      : " << cgc.getPenalty() << endl;
    vector<Rectangle> regions2(2);
    regions2[0] = BBox(0, 0, 0.35, 0.35);
    regions2[1] = BBox(0.45, 0.45, 0.75, 0.75);
    MultRectRegionConstraint cgc1(subPl, regions2);

    cout << "Is satisfied : " << cgc1.isSatisfied() << endl;
    cout << "Penalty      : " << cgc1.getPenalty() << endl;
  }
}
