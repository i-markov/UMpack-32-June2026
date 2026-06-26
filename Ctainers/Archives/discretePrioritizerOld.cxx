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
// discretePrioritizer.cxx: implementation of the DiscretePrioritizer class.
//
//////////////////////////////////////////////////////////////////////
#if (!defined(_MSC_VER) || defined(_INSIDE_DISCRETE_PRIORITIZER_H))

#ifndef _INSIDE_DISCRETE_PRIORITIZER_H
#include "discretePrioritizer.h"
#endif

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////
template <class T>
DiscretePrioritizer<T>::DiscretePrioritizer(vector<const T *> elements,
                                            vector<unsigned> priorities)
    : _eltVec(elements.size(),
              static_cast<DiscretePriorityElement<T> *>(NULL)) {
  abkfatal(elements.size() == priorities.size(), "Mismatched vector sizes");

  // scan for max priority
  _maxPriority = -1;
  unsigned i;
  for (i = 0; i < priorities.size(); i++) {
    if (signed(priorities[i]) > _maxPriority) _maxPriority = priorities[i];
  }

  _heads.clear();
  _tails.clear();
  _heads.insert(_heads.end(), _maxPriority + 1,
                static_cast<DiscretePriorityElement<T> *>(NULL));
  _tails.insert(_tails.end(), _maxPriority + 1,
                static_cast<DiscretePriorityElement<T> *>(NULL));

  for (i = 0; i < elements.size(); i++) {
    unsigned priority = priorities[i];
    const T *pOrigElt = elements[i];
    DiscretePriorityElement<T> *pPrElt =
        new DiscretePriorityElement<T>(i, pOrigElt, priority, _tails[priority]);
    if (!_heads[priority]) _heads[priority] = pPrElt;

    _tails[priority] = pPrElt;

    _eltVec[i] = pPrElt;
  }
}

template <class T>
DiscretePrioritizer<T>::~DiscretePrioritizer() {
  DiscretePrioritizerVecType::iterator iM = _eltVec.begin();
  for (; iM != _eltVec.end(); iM++) {
    delete (*iM);
  }
}

template <class T>
DiscretePriorityElement<T> *DiscretePrioritizer<T>::_disconnect(
    unsigned elementIdx, int newPriority) {
  DiscretePriorityElement<T> *pPrElt = _eltVec[elementIdx];

  if (!pPrElt) return static_cast<DiscretePriorityElement<T> *>(NULL);

  unsigned oldPriority = pPrElt->getPriority();

  DiscretePriorityElement<T> *oldPrev, *oldNext;
  pPrElt->disconnect(oldPrev, oldNext);

  if (!oldPrev) _heads[oldPriority] = oldNext;

  if (!oldNext) _tails[oldPriority] = oldPrev;

  if (newPriority > _maxPriority) {
    _heads.insert(_heads.end(), newPriority - _maxPriority + 1,
                  static_cast<DiscretePriorityElement<T> *>(NULL));
    _tails.insert(_tails.end(), newPriority - _maxPriority + 1,
                  static_cast<DiscretePriorityElement<T> *>(NULL));
  }

  else {  // if we emptied this priority level, scan down for new max
    while (_maxPriority > newPriority && !_heads[_maxPriority]) --_maxPriority;
  }

  return pPrElt;
}

template <class T>
bool DiscretePrioritizer<T>::reQueueToHead(unsigned elementIdx,
                                           int newPriority) {
  DiscretePriorityElement<T> *pPrElt = _disconnect(elementIdx, newPriority);

  if (!pPrElt) return false;  // element not found

  pPrElt->setPriority(newPriority);

  pPrElt->insertBefore(_heads[newPriority]);
  _heads[newPriority] = pPrElt;

  return true;
}
template <class T>
bool DiscretePrioritizer<T>::reQueueToTail(unsigned elementIdx,
                                           int newPriority) {
  DiscretePriorityElement<T> *pPrElt = _disconnect(elementIdx, newPriority);

  if (!pPrElt) return false;  // element not found

  pPrElt->setPriority(newPriority);

  pPrElt->insertAfter(_tails[newPriority]);
  _tails[newPriority] = pPrElt;

  return true;
}

template <class T>
bool DiscretePrioritizer<T>::erase(unsigned elementIdx) {
  DiscretePriorityElement<T> *pPrElt = _disconnect(elementIdx, -1);

  if (!pPrElt) return false;  // element not found

  delete pPrElt;

  _eltVec[elementIdx] = NULL;

  return true;
}

#endif  // !defined(_MSC_VER) || defined(_INSIDE_DISCRETE_PRIORITIZER_H)
