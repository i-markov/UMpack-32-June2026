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

#include "gainBucketHL.h"

using std::endl;
using std::ostream;

/*inline*/ SVGainElement*
BucketArrayHL::getHighest() {  // assumes _highestValid indexes a bucket with
                               // some valid move remaining or -1, and that
                               // lastMoveSuggested is either an element in that
                               // bucket, or NULL if _highestValie == -1
  abkassert(_maxBucketIdx != -1,
            "must call setMaxGain before using the BucketArrayHL");

  return _lastMoveSuggested;
}

// fixed
/*inline*/ void BucketArrayHL::resetAfterGainUpdate() {
  // bring the _highestNonEmpty index down to the first non-empty
  // bucket (it may be too high as a result of incremental gain updates)
  // and clear any invalidations

  abkassert(_maxBucketIdx != -1,
            "must call setMaxGain before using the BucketArrayHL");

  while (_highestNonEmpty != _weightList.begin() &&
         _highestNonEmpty != _weightList.end() &&
         (*this)[*_highestNonEmpty].isEmpty()) {
    abkfatal((*this)[*_highestNonEmpty].isEmpty(),
             "The weight list says false when the bucketArray isnt empty");
    --_highestNonEmpty;
  }
  if (_highestNonEmpty == _weightList.begin() &&
      _highestNonEmpty != _weightList.end() &&
      (*this)[*_highestNonEmpty].isEmpty())
    _highestNonEmpty = _weightList.end();

  _highestValid = _highestNonEmpty;
  if (_highestValid !=
      _weightList.end())  // there is something left in this array
    _lastMoveSuggested = (*this)[*_highestValid]._next;
  else
    _lastMoveSuggested = NULL;
}

// fixed
/*inline*/ bool BucketArrayHL::invalidateBucket() {
  abkassert(_maxBucketIdx != -1,
            "must call setMaxGain before using the BucketArrayHL");

  if (_highestValid ==
      _weightList.end())  // there is nothing left in this array
  {
    _lastMoveSuggested = NULL;
    return true;
  }

  if (_highestValid != _weightList.begin())
    --_highestValid;
  else
    _highestValid = _weightList.end();
  while (_highestValid != _weightList.begin() &&
         _highestValid != _weightList.end() &&
         (*this)[*_highestValid].isEmpty()) {
    --_highestValid;
  }
  if (_highestValid == _weightList.begin() && (*this)[*_highestValid].isEmpty())
    _highestValid = _weightList.end();

  if (_highestValid !=
      _weightList.end())  // there is something left in this array
    _lastMoveSuggested = (*this)[*_highestValid]._next;
  else
    _lastMoveSuggested = NULL;

  return true;
}

// fixed because no need to change it
/*inline*/ bool BucketArrayHL::invalidateElement() {
  abkassert(_maxBucketIdx != -1,
            "must call setMaxGain before using the BucketArrayHL");

  abkassert(_lastMoveSuggested != NULL, "can't invalidate a NULL move");
  _lastMoveSuggested = _lastMoveSuggested->_next;

  if (_lastMoveSuggested == NULL)  // last one in this bucket
    return invalidateBucket();
  else
    return false;
}

// fixed
/*inline*/ void BucketArrayHL::removeAll() {
  _bucketArray.clear();
  _weightList.clear();

  _highestNonEmpty = _weightList.end();
  resetAfterGainUpdate();
}

// fixed
/*inline*/ bool BucketArrayHL::isConsistent() {
  bool consistent = true;
  for (std::unordered_map<int, hash_elementT, std::hash<int>, std::equal_to<int> >::iterator g =
           _bucketArray.begin();
       g != _bucketArray.end(); ++g) {
    if (!g->second.gain_list.isConsistent()) consistent = false;
  }

  return consistent;
}

// NOTE: gain offset is the offset that gives use the bucket
// coresponding to gain ==0.  Note that this is a vector
// of size 2*_gainOffset+1.
// fixed
/*inline*/ void BucketArrayHL::setupForClip() {
  abkassert(isConsistent(), " inconsistent bucket array ");

  (*this)[_gainOffset];  // create the thing
  std::set<int>::iterator g = _weightList.find(_gainOffset);
  abkfatal(g != _weightList.end(),
           "Must have such a weight, it was created above if not");
  for (--g; g != _weightList.begin(); --g) {
    (*this)[_gainOffset].insertBehind((*this)[*g]);
  }
  (*this)[_gainOffset].insertBehind((*this)[*_weightList.begin()]);

  g = _weightList.find(_gainOffset);
  abkfatal(g != _weightList.end(),
           "Must have such a weight, it was created above if not");
  for (++g; g != _weightList.end(); ++g) {
    (*this)[_gainOffset].insertAhead((*this)[*g]);
  }

  SVGainElement* e;
  for (e = (*this)[_gainOffset]._next; e != NULL; e = e->_next) {
    e->_gainOffset = e->_gain;
    e->_gain = 0;
  }
}

// fixed
/*inline*/ ostream& operator<<(ostream& os, const BucketArrayHL& bucket) {
  os << " BucketArrayHL:" << endl;
  os << "  MaxGain:  " << bucket._gainOffset << endl;

  std::set<int>::const_iterator g = bucket._weightList.end();
  for (--g; g != bucket._weightList.begin(); --g) {
    int weightIdx = *g;
    os << weightIdx - bucket._gainOffset << ") " << (bucket[weightIdx]);
  }
  os << (*bucket._weightList.begin()) - bucket._gainOffset << ") "
     << bucket[(*bucket._weightList.begin())];

  os << endl;
  return os;
}

// fixed
/*inline*/ void BucketArrayHL::setMaxGain(unsigned maxG) {
  _gainOffset = maxG;
  _maxBucketIdx = 2 * _gainOffset;

  _highestNonEmpty = _highestValid = _weightList.end();
}

// fixed
/*inline*/ void BucketArrayHL::push(SVGainElement& elmt) {
  abkassert(_maxBucketIdx != -1,
            "must call setMaxGain before using the BucketArrayHL");

  int bucket = std::max(std::min(_maxBucketIdx, elmt._gain + _gainOffset), 0);
  (*this)[bucket].push(elmt);
  if (_highestNonEmpty == _weightList.end() || *_highestNonEmpty < bucket)
    _highestNonEmpty = _weightList.find(bucket);
}
