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

#ifndef __PARTEVAL_H_
#define __PARTEVAL_H_

// interfaces
#include "evalRegistry.h"
#include "partEvalXFace.h"
#include "univPartEval.h"

// reusable implementations of generic features
#include "netTallies.h"
#include "netVec.h"
#include "talliesWConfigIds.h"
#include "talliesWCosts.h"

// specific evaluators
#include "bbox1dim.h"
#include "netCutWBits.h"
#include "netCutWConfigIds.h"
#include "netCutWNetVec.h"
/*
#include "bbox1dimWCheng.h"
#include "bbox2dim.h"
#include "bbox2dimWCheng.h"
#include "bbox2dimWRSMT.h"
#include "hbbox.h"
#include "hbbox0.h"
#include "hbbox0wCheng.h"
#include "hbbox0wRSMT.h"
#include "hbboxWCheng.h"
#include "hbboxWRSMT.h"
*/
#include "netCut2way.h"
#include "netCut2wayWWeights.h"
#include "strayNodes.h"
#include "strayNodes2way.h"

#endif
