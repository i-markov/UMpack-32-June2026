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
// prioritizer.cxx: implementation of the Prioritizer class.
//
//////////////////////////////////////////////////////////////////////
#if (!defined(_MSC_VER) || defined(_INSIDE_PRIORITIZER_H))

#ifndef _INSIDE_PRIORITIZER_H
#include "prioritizer.h"
#endif

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////
template <class T>
Prioritizer<T>::Prioritizer(vector<const T *> elements,
                            vector<unsigned> priorities, bool populateMap)
    : _mapPopulated(populateMap),
      _eltVec(elements.size(), static_cast<PriorityElement<T> *>(NULL)) {
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
                static_cast<PriorityElement<T> *>(NULL));
  _tails.insert(_tails.end(), _maxPriority + 1,
                static_cast<PriorityElement<T> *>(NULL));

  for (i = 0; i < elements.size(); i++) {
    unsigned priority = priorities[i];
    const T *pOrigElt = elements[i];
    PriorityElement<T> *pPrElt =
        new PriorityElement<T>(i, pOrigElt, priority, _tails[priority]);
    if (!_heads[priority]) _heads[priority] = pPrElt;

    _tails[priority] = pPrElt;

    if (_mapPopulated) {
      pair<T const *const, unsigned> tempPair(pOrigElt, i);
      pair<MikePrioMapType::iterator, bool> success = _eltMap.insert(tempPair);

      abkfatal(success.second, "Duplicate objects passed to Prioritizer");
    }

    _eltVec[i] = pPrElt;
  }
}

template <class T>
Prioritizer<T>::~Prioritizer() {
  MikePrioVecType::iterator iM = _eltVec.begin();
  for (; iM != _eltVec.end(); iM++) {
    delete (*iM);
  }
}

template <class T>
unsigned Prioritizer<T>::getEltIndex(const T *element) const {
  abkfatal(_mapPopulated,
           "Attempt to find element from unpopulated "
           "map");
  MikePrioMapType::const_iterator iM;
  iM = _eltMap.find(element);

  if (iM == _eltMap.end()) return UINT_MAX;  // element not found

  return (*iM).second;
}

template <class T>
PriorityElement<T> *Prioritizer<T>::_disconnect(unsigned elementIdx,
                                                int newPriority) {
  PriorityElement<T> *pPrElt = _eltVec[elementIdx];

  if (!pPrElt) return static_cast<PriorityElement<T> *>(NULL);

  unsigned oldPriority = pPrElt->getPriority();

  PriorityElement<T> *oldPrev, *oldNext;
  pPrElt->disconnect(oldPrev, oldNext);

  if (!oldPrev) _heads[oldPriority] = oldNext;

  if (!oldNext) _tails[oldPriority] = oldPrev;

  if (newPriority > _maxPriority) {
    _heads.insert(_heads.end(), newPriority - _maxPriority + 1,
                  static_cast<PriorityElement<T> *>(NULL));
    _tails.insert(_tails.end(), newPriority - _maxPriority + 1,
                  static_cast<PriorityElement<T> *>(NULL));
  }

  else {  // if we emptied this priority level, scan down for new max
    while (_maxPriority > newPriority && !_heads[_maxPriority]) --_maxPriority;
  }

  return pPrElt;
}

template <class T>
bool Prioritizer<T>::reQueueToHead(const T *element, int newPriority) {
  unsigned elementIdx = getEltIndex(element);

  if (elementIdx == UINT_MAX) return false;

  return reQueueToHead(elementIdx, newPriority);
}
template <class T>
bool Prioritizer<T>::reQueueToTail(const T *element, int newPriority) {
  unsigned elementIdx = getEltIndex(element);

  if (elementIdx == UINT_MAX) return false;

  return reQueueToTail(elementIdx, newPriority);
}

template <class T>
bool Prioritizer<T>::erase(const T *element) {
  unsigned elementIdx = getEltIndex(element);

  if (elementIdx == UINT_MAX) return false;

  bool retval = erase(elementIdx);

  _eltMap.erase(element);

  return retval;
}
template <class T>
bool Prioritizer<T>::reQueueToHead(unsigned elementIdx, int newPriority) {
  PriorityElement<T> *pPrElt = _disconnect(elementIdx, newPriority);

  if (!pPrElt) return false;  // element not found

  pPrElt->setPriority(newPriority);

  pPrElt->insertBefore(_heads[newPriority]);
  _heads[newPriority] = pPrElt;

  return true;
}
template <class T>
bool Prioritizer<T>::reQueueToTail(unsigned elementIdx, int newPriority) {
  PriorityElement<T> *pPrElt = _disconnect(elementIdx, newPriority);

  if (!pPrElt) return false;  // element not found

  pPrElt->setPriority(newPriority);

  pPrElt->insertAfter(_tails[newPriority]);
  _tails[newPriority] = pPrElt;

  return true;
}

template <class T>
bool Prioritizer<T>::erase(unsigned elementIdx) {
  PriorityElement<T> *pPrElt = _disconnect(elementIdx, -1);

  if (!pPrElt) return false;  // element not found

  delete pPrElt;

  _eltVec[elementIdx] = NULL;

  return true;
}

#endif  // !defined(_MSC_VER) || defined(_INSIDE_PRIORITIZER_H)
