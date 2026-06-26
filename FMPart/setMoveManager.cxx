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

#include "dcGainCont.h"
#include "fmMoveManagerNC2w.h"
#include "fmNetCut2WayMM.h"
#include "fmPart.h"

void FMPartitioner::setMoveManagerAndSwitchBox() {
  unsigned numPart = _problem.getNumPartitions();
  abkfatal2(numPart >= 2, " Too few partitions : ", numPart);

  _activeMoveMgr = _moveMgrLIFO;

  if (!_params.useClip &&  // D&D lifo is ONLY D&D. otherwise, it uses both
      MoveManagerType::Type(_params.moveMan) == MoveManagerType::FMDD)
    _switchBox = new mmSwitchBoxNoSwaps();
  else if (_params.doFirstLIFOPass)
    _switchBox = new mmSwitchBoxConservative();
  else
    _switchBox = new mmSwitchBoxSwapOnce();

  switch (_params.moveMan) {
    case MoveManagerType::FM:
      _moveMgrLIFO = new FMNetCut2WayMoveManager(_problem, _params,
                                                 _params.allowCorkingNodes);
      _moveMgr2 = NULL;

      break;
    case MoveManagerType::FMDD:
      _moveMgrLIFO = new FMNetCut2WayMoveManager(_problem, _params,
                                                 _params.allowCorkingNodes);
      _moveMgr2 = _moveMgrLIFO;
    default:
      abkfatal(1, "This piece of code should not be reached");
  };

  abkfatal(_moveMgrLIFO, " Internal error: LIFO move manager not initialized");

  _moveMgrLIFO->randomized = _params.randomized;
  _moveMgrLIFO->unCorking = _params.unCorking;
  _moveMgrLIFO->allowIllegalMoves = _params.allowIllegalMoves;

  if (_moveMgr2) {
    _moveMgr2->unCorking = _params.unCorking;
    _moveMgr2->randomized = _params.randomized;
    _moveMgr2->allowIllegalMoves = _params.allowIllegalMoves;
    _moveMgr2->wiggleTerms = _params.wiggleTerms;
  }
}

void FMPartitioner::swapMoveManagers() {
  delete _activeMoveMgr;
  if (_moveMgrLIFO == NULL) {
    _moveMgr2 = NULL;
    _moveMgrLIFO = new FMNetCut2WayMoveManager(_problem, _params,
                                               _params.allowCorkingNodes);

    _moveMgrLIFO->randomized = _params.randomized;
    _moveMgrLIFO->unCorking = _params.unCorking;
    _moveMgrLIFO->allowIllegalMoves = _params.allowIllegalMoves;

    _activeMoveMgr = _moveMgrLIFO;
  } else {
    _moveMgrLIFO = NULL;
    _moveMgr2 =
        new FMMoveManagerNC2w(_problem, _params, _params.allowCorkingNodes);

    _moveMgr2->unCorking = _params.unCorking;
    _moveMgr2->randomized = _params.randomized;
    _moveMgr2->allowIllegalMoves = _params.allowIllegalMoves;
    _moveMgr2->wiggleTerms = _params.wiggleTerms;

    _activeMoveMgr = _moveMgr2;
  }
}
