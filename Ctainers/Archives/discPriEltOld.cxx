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
// discPriElt.cxx: implementation of the DiscretePriorityElement class.
//
//////////////////////////////////////////////////////////////////////
#if (!defined(_MSC_VER) || defined(_INSIDE_DISCPRIELT_H))

#ifndef _INSIDE_DISCPRIELT_H
#include "discPriElt.h"
#endif  //_INSIDE_DISCPRIELT_H

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

template <class T>
DiscretePriorityElement<T>::DiscretePriorityElement(
    unsigned index, const T *bareElement, int priority,
    DiscretePriorityElement *prevElt)
    : _index(index),
      _bareElement(bareElement),
      _prev(prevElt),
      _next(NULL),
      _priority(priority) {
  if (prevElt) prevElt->_next = this;
}

template <class T>
DiscretePriorityElement<T>::~DiscretePriorityElement() {}

template <class T>
void DiscretePriorityElement<T>::disconnect(DiscretePriorityElement *&oldPrev,
                                            DiscretePriorityElement *&oldNext) {
  oldPrev = _prev;
  oldNext = _next;
  if (_prev) _prev->_next = _next;
  if (_next) _next->_prev = _prev;
  _prev = _next = NULL;
}

template <class T>
void DiscretePriorityElement<T>::insertBefore(
    DiscretePriorityElement *nextElt) {
  abkfatal(!_next && !_prev,
           "attempt to insert already-connected "
           "element");

  _next = nextElt;

  if (nextElt) {
    _prev = nextElt->_prev;
    nextElt->_prev = this;
  }

  if (_prev) _prev->_next = this;
}

template <class T>
void DiscretePriorityElement<T>::insertAfter(DiscretePriorityElement *prevElt) {
  abkfatal(!_next && !_prev,
           "attempt to insert already-connected "
           "element");

  _prev = prevElt;

  if (prevElt) {
    _next = prevElt->_next;
    prevElt->_next = this;
  }

  if (_next) _next->_prev = this;
}

#endif  // !defined(_MSC_VER) || defined(_INSIDE_DISCPRIELT_H)
