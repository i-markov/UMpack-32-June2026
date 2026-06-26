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
//! author="Mike Oliver 1999.4.27"
//! CONTACTS="abk@cs.ucsd.edu, igor.markov1@gmail.com"
// prioritizer.h: interface for the Prioritizer class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(_PRIORITIZER_H_INCLUDED_)
#define _PRIORITIZER_H_INCLUDED_

#include <map.h>
#include <vector.h>

#define MikePrioMapType map<const T *, unsigned, less<const T *> >
#define MikePrioVecType vector<PriorityElement<T> *>

#if _MSC_VER > 1000
#pragma once
#endif  // _MSC_VER > 1000

#include "priElt.h"

template <class T>
class Prioritizer {
 protected:
  MikePrioMapType _eltMap;
  MikePrioVecType _eltVec;

  // The following arrays are parallel.  _heads[i] points to
  // the beginning of the linked list with priority i; _tails[i]
  // points to the last element of that list.
  vector<PriorityElement<T> *> _heads;
  vector<PriorityElement<T> *> _tails;
  int _maxPriority;

  PriorityElement<T> *_disconnect(unsigned elementIdx, int newPriority);

 public:
  Prioritizer(vector<const T *> elements, vector<unsigned> priorities,
              bool populateMap = false);
  virtual ~Prioritizer();

  const T *top() const {
    if (_maxPriority < 0)
      return NULL;
    else
      return _heads[_maxPriority]->getBareElement();
  }

  unsigned topEltIndex() const {
    if (_maxPriority < 0)
      return UINT_MAX;
    else
      return _heads[_maxPriority]->getIndex();
  }

  bool reQueueToHead(unsigned elementIdx, int newPriority);
  bool reQueueToHead(const T *element, int newPriority);
  bool reQueueToTail(unsigned elementIdx, int newPriority);
  bool reQueueToTail(const T *element, int newPriority);

  bool erase(unsigned elementIdx);
  bool erase(const T *element);

  unsigned getEltIndex(const T *element) const;
};

#ifdef _MSC_VER
#define _INSIDE_PRIORITIZER_H
#include "prioritizer.cxx"
#undef _INSIDE_PRIORITIZER_H
#endif

#endif  // !defined(_PRIORITIZER_H_INCLUDED_)
