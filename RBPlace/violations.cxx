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

// 990107  aec - complete restructuring of class implementation and interface

#ifdef _MSC_VER
#pragma warning(disable : 4786)
#endif

#include "violations.h"

#include <cstring>

#include "ABKCommon/abkassert.h"

using std::ostream;

Violations::Violations(const char* txt) : _violations(0) {
  if (strstr(txt, "offSite"))
    _violations += OffSite;
  else if (strstr(txt, "hOverlap"))
    _violations += HOverlap;
  else if (strstr(txt, "vOverlap"))
    _violations += VOverlap;
  else if (strstr(txt, "wOverlap"))
    _violations += WOverlap;
  else if (strstr(txt, "rOrient"))
    _violations += ROrient;
  else if (strstr(txt, "cOrient"))
    _violations += COrient;
  else if (strstr(txt, "subRow"))
    _violations += SubRow;
  else if (strstr(txt, "offRow"))
    _violations += OffRow;
  else
    abkfatal2(0, "invalid string argument to Violations ctor:", txt);
}

ostream& operator<<(ostream& os, const Violations& v) {
  if (!v.hasViolations())
    os << "<NONE>";
  else {
    os << "<";
    if (v.hasViolation(Violations::OffSite)) os << " offsite";
    if (v.hasViolation(Violations::HOverlap)) os << " hOverlap";
    if (v.hasViolation(Violations::VOverlap)) os << " vOverlap";
    if (v.hasViolation(Violations::WOverlap)) os << " wOVerlap";
    if (v.hasViolation(Violations::ROrient)) os << " rOrient";
    if (v.hasViolation(Violations::COrient)) os << " cOrient";
    if (v.hasViolation(Violations::SubRow)) os << " subRow";
    if (v.hasViolation(Violations::OffRow)) os << " offRow";
    os << ">";
  }

  return os;
}
