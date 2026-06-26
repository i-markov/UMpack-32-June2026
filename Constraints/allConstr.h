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

//! author =" Aug 16, 1997 by Igor Markov"

#ifndef _ALLCONSTRS_H_
#define _ALLCONSTRS_H_

/* if you don't need all of these headers, better include them individually
 * (see also dbConstr.h) */

#include "fxTypeConstr.h"
#include "otherConstr.h"
#include "regionConstr.h"

/* descriptions *follow* class declarations */

/*  FIXED TYPE CONSTRAINTS  */
class FixedGroupConstraint;
// all cells must have x_coord and y_coord fixed at given points
class FixedXGroupConstraint;
// all cells must have x_coord fixed at given values
class FixedYGroupConstraint;
// all cells must have y_coord fixed at given values

/*  RECTILINEAR REGION CONSTRAINTS */

class SoftGroupConstraint;
// the dimensions of the BBox of given cells are bounded
// (see LEF/DEF manual p.6-51.)
class RectRegionConstraint;
// all ref'ed cells must be located in the rectangular region
class MultRectRegionConstraint;
// all ref'ed cells must be located in of of given rectangular regions
class UtilRectRegionConstraint;
// The rectangular region must satisfy utilization constraints
/*  OTHER CONSTRAINTS */

class TetheredGroupConstraint;
// all ref'ed cells must be tethered (to same tolerance dist)
//  at their locations
class COGConstraint;
// referenced cells must have COG at the given center point
class EqualXConstraint;
// referenced cells must have equal x coordinates
class EqualYConstraint;
// referenced cells must have equal y coordinates
class StayTogetherConstraint;
// referenced cells must be located in a given
// rectangular region centered at their center of gravity

#endif
