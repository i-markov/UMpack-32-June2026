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

// Created: July 25, 1997 by Igor Markov
// 970819 mro parenthesized references in _pullBack()
// 970828 ilm separated restriction from pull-back
// 970828 ilm reimplemented restriction with partial inverse
//                   (was with sorting and binary_search)

#include <ABKCommon/abkassert.h>
#include <Constraints/constraint.h>
#include <Placement/placement.h>
#include <Placement/subPlacement.h>

#include <algorithm>
#include <iostream>
#include <vector>

using std::cerr;
using std::vector;

const char* Constraint::Type::_text = NULL;

Constraint::Type::operator const char*() const {
  switch (_type) {
    case FixedGroup:
      _text = "Fixed   Group Constraint";
      break;
    case FixedXGroup:
      _text = "Fixed X Group Constraint";
      break;
    case FixedYGroup:
      _text = "Fixed Y Group Constraint";
      break;
    case SoftGroup:
      _text = "  Soft  Group Constraint";
      break;
    case RectRegion:
      _text = " Rect. Region Constraint";
      break;
    case MultRectRegion:
      _text = "Mult Rect. Region Constraint";
      break;
    case COG:
      _text = "Ctr of Grvity Constraint";
      break;
    case TetheredGroup:
      _text = " Tethered Grp Constraint";
      break;
    case EqualX:
      _text = "      Equal X Constraint";
      break;
    case EqualY:
      _text = "      Equal Y Constraint";
      break;
    case StayTogether:
      _text = "Stay Together Constraint";
      break;
    case UtilRectRegion:
      _text = "Utilization Rect. Region Constraint";
      break;
    default:
      cerr << "Unknown constraint Type " << _type << " ";
      abkfatal(0, "");
  }
  return _text;
}

/* Old restriction with sorting and binary_search : obsolete
Constraint* Constraint::_restriction(const vector<unsigned>& sortedIdx) const
{
    const vector<unsigned> (&splIndices) (_spl.getMapping());
    vector<unsigned> hitIndices, hitRealIndices;

    for(unsigned k=0;k<splIndices.size();k++)
      if (binary_search(sortedIdx.begin(),sortedIdx.end(),splIndices[k]))
      {
         hitIndices.push_back(k);
         hitRealIndices.push_back(splIndices[k]);
      }

    if (hitIndices.size()==0) return NULL;
    else
    {
        Mapping mp(hitIndices.size(), _spl.getPlacement().getSize(),
                   hitRealIndices);
        SubPlacement newSubPl(_spl.getPlacement(),mp);
        return _virtPullBack(hitIndices,newSubPl);
    }
}
*/

// New restriction with partial inverse
Constraint* Constraint::_restriction(
    const std::vector<unsigned>& partialInv) const {
  const std::vector<unsigned>& splIndices = _spl.getMapping().getVector();
  vector<unsigned> hitIndices;
  std::vector<unsigned> hitRealIndices;

  for (unsigned k = 0; k < splIndices.size(); k++)
    if (partialInv[splIndices[k]] != unsigned(-1)) {
      hitIndices.push_back(k);
      hitRealIndices.push_back(splIndices[k]);
    }

  if (hitIndices.size() == 0)
    return NULL;
  else {
    Mapping mp(hitIndices.size(), _spl.getPlacement().getSize(),
               hitRealIndices);
    SubPlacement newSubPl(_spl.getPlacement(), mp);
    return _virtPullBack(hitIndices, newSubPl);
  }
}

Constraint* Constraint::_pullBack(const std::vector<unsigned>& partalInverse,
                                  const Placement& pullBackPl) const {
  const std::vector<unsigned>& splIndices = _spl.getMapping().getVector();
  vector<unsigned> hitIndices;
  std::vector<unsigned> hitRealIndices;

  for (unsigned k = 0; k < splIndices.size(); k++)
    if (partalInverse[splIndices[k]] != unsigned(-1)) {
      hitIndices.push_back(k);
      hitRealIndices.push_back(partalInverse[splIndices[k]]);
    }

  if (hitIndices.size() == 0)
    return NULL;
  else {
    Mapping mp(hitIndices.size(), pullBackPl.getSize(), hitRealIndices);
    SubPlacement newSubPl(pullBackPl, mp);
    return _virtPullBack(hitIndices, newSubPl);
  }
}
Constraint::~Constraint() {}
