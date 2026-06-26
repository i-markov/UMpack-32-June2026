/**************************************************************************
***    
*** Copyright (c) 2026 IonQ, Inc.
*** Copyright (c) 1995-2000 Regents of the University of California,
***               Andrew E. Caldwell, Andrew B. Kahng and Igor L. Markov
*** Copyright (c) 2000-2012 Regents of the University of Michigan,
***               Saurabh N. Adya, Jarrod A. Roy, David A. Papa and
***               Igor L. Markov
***
***  Contact author(s): abk@cs.ucsd.edu, imarkov@umich.edu
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


// Created 990326 by Stefanus Mantik, VLSI CAD ABKGROUP, UCLA/CS

#ifndef _HEAP_GAIN_CONTAINER_H
#define _HEAP_GAIN_CONTAINER_H

#include <algorithm>
#include "Ctainers/bitBoard.h"

class HeapGainCtrLessComp
{
    std::vector<double>&          _gain;
  public:
    HeapGainCtrLessComp(std::vector<double>& gain) : _gain(gain) {}
    bool operator()(unsigned x, unsigned y)
    { return _gain[x] < _gain[y]; }
};

class HeapGainContainer
{
    std::vector<double>&   _gain;
    std::vector<unsigned>  _idxToHeap;
    std::vector<unsigned>  _heap;

  public:

    HeapGainContainer(std::vector<double>& gain, const BitBoard& bb) :
        _gain(gain), _idxToHeap(gain.size(), UINT_MAX),
        _heap(bb.getNumBitsSet(), UINT_MAX)
    {
        const std::vector<unsigned>& indices = bb.getIndicesOfSetBits();
        unsigned i;
        for(i = 0; i < indices.size(); i++)
            _heap[i] = indices[i];
        std::make_heap(_heap.begin(), _heap.end(), HeapGainCtrLessComp(_gain));
        for(i = 0; i < _heap.size(); i++)
            _idxToHeap[_heap[i]] = i;
    }

    void heapify(unsigned index)
    {
        abkfatal(index < _gain.size(), "index out of range");
        bool moved = false;
        if(_heap.size() <= 1) return;
        if(_idxToHeap[index] > 0)    // upward heapify
        {
            unsigned currentId = _idxToHeap[index];
            while(currentId > 0)
            {
                unsigned parentId = 
			static_cast<unsigned>(floor((currentId - 1)/2.0));
                if(_gain[_heap[parentId]] < _gain[_heap[currentId]])
                {
                    //swap(_heap[parentId], _heap[currentId]);
                    unsigned swap_temp = _heap[parentId];
                    _heap[parentId] = _heap[currentId];
                    _heap[currentId] = swap_temp;
                    _idxToHeap[_heap[parentId]] = parentId;
                    _idxToHeap[_heap[currentId]] = currentId;
                    currentId = parentId;
                    moved = true;
                }
                else break;
            }
        }
        if(moved) return;
        if(_idxToHeap[index] < _heap.size() - 1)    // downward heapify
        {
            unsigned currentId = _idxToHeap[index];
            while(true)
            {
                unsigned leftId = (currentId * 2) + 1,
                         rightId = (currentId * 2) + 2;
                if(leftId >= _heap.size()) break;
                unsigned nextId;
                if(rightId >= _heap.size()) nextId = leftId;
                else if(_gain[_heap[rightId]] < _gain[_heap[leftId]])
                    nextId = leftId;
                else
                    nextId = rightId;
                if(_gain[_heap[currentId]] < _gain[_heap[nextId]])
                {
                    //swap(_heap[currentId], _heap[nextId]);
                    unsigned swap_temp = _heap[currentId];
                    _heap[currentId] = _heap[nextId];
                    _heap[nextId] = swap_temp;
                    _idxToHeap[_heap[nextId]] = nextId;
                    _idxToHeap[_heap[currentId]] = currentId;
                    currentId = nextId;
                }
                else break;
            }
        }
    }

    unsigned getTop() { return _heap[0]; }
};

#endif
