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
//
//////////////////////////////////////////////////////////////////////
// discPriElt.h: interface for the DiscretePriorityElement class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(_DISCPRIELT_H__INCLUDED_)
#define _DISCPRIELT_H__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif  // _MSC_VER > 1000

#include <abkcommon.h>

template <class T>
class DiscretePriorityElement {
 protected:
  const unsigned _index;
  const T *_bareElement;
  DiscretePriorityElement *_prev;
  DiscretePriorityElement *_next;
  int _priority;

 public:
  DiscretePriorityElement(unsigned index, const T *bareElement, int priority,
                          DiscretePriorityElement *prevElt = NULL);
  virtual ~DiscretePriorityElement();

  unsigned getPriority() const { return _priority; }
  void setPriority(int pr) { _priority = pr; }

  unsigned getIndex() const { return _index; }

  const T *getBareElement() const { return _bareElement; }

  void disconnect(DiscretePriorityElement *&oldPrev,
                  DiscretePriorityElement *&oldNext);

  void insertBefore(DiscretePriorityElement *nextElt);
  void insertAfter(DiscretePriorityElement *prevElt);
};
#ifdef _MSC_VER
#define _INSIDE_DISCPRIELT_H
#include "discPriElt.cxx"
#undef _INSIDE_DISCPRIELT_H
#endif

#endif  // !defined(_DISCPRIELT_H__INCLUDED_)
