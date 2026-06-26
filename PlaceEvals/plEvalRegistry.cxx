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

// created by Stefanus Mantik on 06/19/98
#ifdef _MSC_VER
#pragma warning(disable : 4786)
#endif

#include "plEvalRegistry.h"

#include <limits>
#include <strings.h>

#include "ABKCommon/paramproc.h"
#include "HGraph/hgFixed.h"
#include "Placement/layoutBBoxes.h"
#include "Placement/placement.h"
#include "placeEvalXFace.h"

// Now include headers for all evaluators available
#include "bboxPlEval.h"
#include "bboxPlEvalWt.h"
#include "bboxWChengPlEval.h"
#include "bboxWRSMTPlEval.h"
#include "hbbox0PlEval.h"
#include "hbbox0wChengPlEval.h"
#include "hbbox0wRSMTPlEval.h"
#include "hbbox6PlEval.h"
#include "hbbox6wRSMTPlEval.h"
#include "hbboxNPlEval.h"
#include "hbboxNwChengPlEval.h"
#include "hbboxNwRSMTPlEval.h"
#include "hbboxPlEval.h"
#include "hbboxWChengPlEval.h"
#include "hbboxWRSMTPlEval.h"

using std::cout;
using std::endl;
using std::ostream;
using std::vector;

#define SEARCH_FOR_PLEVAL(NAME) \
  if (!strcasecmp(evalName, #NAME)) _type = NAME;

PlaceEvalType::PlaceEvalType(int argc, const char* argv[]) {
  _type = static_cast<Type>(std::numeric_limits<unsigned>::max());
  StringParam evalName("placeEval", argc, argv);
  StringParam helpRequest1("h", argc, argv);
  StringParam helpRequest2("help", argc, argv);

  if (helpRequest1.found() || helpRequest2.found()) {
    cout << " -placeEvals  BBoxPlEval  | BBoxPlEvalWt | BBoxWCheng | BBoxWRSMT "
         << "          | HBBox | HBBoxWCheng | HBBoxWRSMT "
         << "          | HBBox0 | HBBox0wCheng | HBBox0wRSMT "
         << "          | HBBoxN | HBBoxNwCheng | HBBoxNwRSMT"
         << "          | HBBox6 | HBBox6wRSMT" << endl;
    return;
  }

  if (evalName.found()) {
    SEARCH_FOR_PLEVAL(BBoxPlEval)
    else SEARCH_FOR_PLEVAL(BBoxPlEvalWt) else SEARCH_FOR_PLEVAL(BBoxWCheng) else SEARCH_FOR_PLEVAL(BBoxWRSMT) else SEARCH_FOR_PLEVAL(
        HBBox) else SEARCH_FOR_PLEVAL(HBBoxWCheng) else SEARCH_FOR_PLEVAL(HBBoxWRSMT) else SEARCH_FOR_PLEVAL(HBBox0) else SEARCH_FOR_PLEVAL(HBBox0wCheng) else SEARCH_FOR_PLEVAL(HBBox0wRSMT) else SEARCH_FOR_PLEVAL(HBBoxN) else SEARCH_FOR_PLEVAL(HBBoxNwCheng) else SEARCH_FOR_PLEVAL(HBBoxNwRSMT) else SEARCH_FOR_PLEVAL(HBBox6) else SEARCH_FOR_PLEVAL(HBBox6wRSMT) else abkfatal3(0,
                                                                                                                                                                                                                                                                                                                                                                                        " -placeEval followed by an unkown evaluator name \"",
                                                                                                                                                                                                                                                                                                                                                                                        evalName,
                                                                                                                                                                                                                                                                                                                                                                                        "\" in command line\n");
  } else {
    cout << "No -placeEval <name> in command line; assuming BBoxPlEval ... "
         << endl;
    _type = BBoxPlEval;
  }
}

#define COMPUTE_COST_PLEVAL(hgraph, place, EVAL) \
  {                                              \
    ::EVAL eval(hgraph, place);                  \
    return eval.getCost();                       \
  }

double PlaceEvalType ::computeCost(const HGraphFixed& hgraph,
                                   const Placement& place) {
  switch (_type) {
    case BBoxPlEval:
      COMPUTE_COST_PLEVAL(hgraph, place, BBoxPlEval)
    case BBoxPlEvalWt:
      COMPUTE_COST_PLEVAL(hgraph, place, BBoxPlEvalWt)
    case BBoxWCheng:
      COMPUTE_COST_PLEVAL(hgraph, place, BBoxWChengPlEval)
    case BBoxWRSMT:
      COMPUTE_COST_PLEVAL(hgraph, place, BBoxWRSMTPlEval)
    case HBBox:
    case HBBoxWCheng:
    case HBBoxWRSMT:
    case HBBox0:
    case HBBox0wCheng:
    case HBBox0wRSMT:
    case HBBoxN:
    case HBBoxNwCheng:
    case HBBoxNwRSMT:
    case HBBox6:
    case HBBox6wRSMT:
      abkfatal(0, "Require regions list and assignment.");
    default:
      abkfatal(0, " Evaluator requested unknown to PlaceEvalType class\n");
  }
  return UINT_MAX;
}

#define COMPUTE_COST_WREGS_PLEVAL(hgraph, regions, assignment, EVAL) \
  {                                                                  \
    ::EVAL eval(hgraph, regions, assignment);                        \
    return eval.getCost();                                           \
  }

double PlaceEvalType ::computeCost(const HGraphFixed& hgraph,
                                   const LayoutBBoxes& regions,
                                   vector<unsigned>& assignment) {
  switch (_type) {
    case BBoxPlEval:
    case BBoxPlEvalWt:
    case BBoxWCheng:
    case BBoxWRSMT:
      abkfatal(0, "Require placement.");
    case HBBox:
      COMPUTE_COST_WREGS_PLEVAL(hgraph, regions, assignment, HBBoxPlEval)
    case HBBoxWCheng:
      COMPUTE_COST_WREGS_PLEVAL(hgraph, regions, assignment, HBBoxWChengPlEval)
    case HBBoxWRSMT:
      COMPUTE_COST_WREGS_PLEVAL(hgraph, regions, assignment, HBBoxWRSMTPlEval)
    case HBBox0:
      COMPUTE_COST_WREGS_PLEVAL(hgraph, regions, assignment, HBBox0PlEval)
    case HBBox0wCheng:
      COMPUTE_COST_WREGS_PLEVAL(hgraph, regions, assignment, HBBox0wChengPlEval)
    case HBBox0wRSMT:
      COMPUTE_COST_WREGS_PLEVAL(hgraph, regions, assignment, HBBox0wRSMTPlEval)
    case HBBoxN:
      COMPUTE_COST_WREGS_PLEVAL(hgraph, regions, assignment, HBBoxNPlEval)
    case HBBoxNwCheng:
      COMPUTE_COST_WREGS_PLEVAL(hgraph, regions, assignment, HBBoxNwChengPlEval)
    case HBBoxNwRSMT:
      COMPUTE_COST_WREGS_PLEVAL(hgraph, regions, assignment, HBBoxNwRSMTPlEval)
    case HBBox6:
      COMPUTE_COST_WREGS_PLEVAL(hgraph, regions, assignment, HBBox6PlEval)
    case HBBox6wRSMT:
      COMPUTE_COST_WREGS_PLEVAL(hgraph, regions, assignment, HBBox6wRSMTPlEval)
    default:
      abkfatal(0, " Evaluator requested unknown to PlaceEvalType class\n");
  }
  return UINT_MAX;
}

#define COMPUTE_COST_WRPL_PLEVAL(hgraph, regions, placement, EVAL) \
  {                                                                \
    ::EVAL eval(hgraph, regions, placement);                       \
    return eval.getCost();                                         \
  }

double PlaceEvalType ::computeCost(const HGraphFixed& hgraph,
                                   const LayoutBBoxes& regions,
                                   const Placement& placement) {
  switch (_type) {
    // case BBoxPlEval :  case BBoxWCheng :  case BBoxWRSMT  :
    // abkfatal(0, "Require placement.");
    case BBoxPlEval:
      COMPUTE_COST_WRPL_PLEVAL(hgraph, regions, placement, BBoxPlEval)
    case BBoxPlEvalWt:
      COMPUTE_COST_WRPL_PLEVAL(hgraph, regions, placement, BBoxPlEvalWt)
    case BBoxWCheng:
      COMPUTE_COST_WRPL_PLEVAL(hgraph, regions, placement, BBoxWChengPlEval)
    case BBoxWRSMT:
      COMPUTE_COST_WRPL_PLEVAL(hgraph, regions, placement, BBoxWRSMTPlEval)
    case HBBox:
      COMPUTE_COST_WRPL_PLEVAL(hgraph, regions, placement, HBBoxPlEval)
    case HBBoxWCheng:
      COMPUTE_COST_WRPL_PLEVAL(hgraph, regions, placement, HBBoxWChengPlEval)
    case HBBoxWRSMT:
      COMPUTE_COST_WRPL_PLEVAL(hgraph, regions, placement, HBBoxWRSMTPlEval)
    case HBBox0:
      COMPUTE_COST_WRPL_PLEVAL(hgraph, regions, placement, HBBox0PlEval)
    case HBBox0wCheng:
      COMPUTE_COST_WRPL_PLEVAL(hgraph, regions, placement, HBBox0wChengPlEval)
    case HBBox0wRSMT:
      COMPUTE_COST_WRPL_PLEVAL(hgraph, regions, placement, HBBox0wRSMTPlEval)
    case HBBoxN:
      COMPUTE_COST_WRPL_PLEVAL(hgraph, regions, placement, HBBoxNPlEval)
    case HBBoxNwCheng:
      COMPUTE_COST_WRPL_PLEVAL(hgraph, regions, placement, HBBoxNwChengPlEval)
    case HBBoxNwRSMT:
      COMPUTE_COST_WRPL_PLEVAL(hgraph, regions, placement, HBBoxNwRSMTPlEval)
    case HBBox6:
      COMPUTE_COST_WRPL_PLEVAL(hgraph, regions, placement, HBBox6PlEval)
    case HBBox6wRSMT:
      COMPUTE_COST_WRPL_PLEVAL(hgraph, regions, placement, HBBox6wRSMTPlEval)
    default:
      abkfatal(0, " Evaluator requested unknown to PlaceEvalType class\n");
  }
  return UINT_MAX;
}

#define PRINT_COST_PLEVAL(hgraph, place, EVAL) \
  {                                            \
    ::EVAL eval(hgraph, place);                \
    cout << eval;                              \
    break;                                     \
  }

void PlaceEvalType ::costPerNet(const HGraphFixed& hgraph,
                                const Placement& place) {
  switch (_type) {
    case BBoxPlEval:
      PRINT_COST_PLEVAL(hgraph, place, BBoxPlEval)
    case BBoxPlEvalWt:
      PRINT_COST_PLEVAL(hgraph, place, BBoxPlEvalWt)
    case BBoxWCheng:
      PRINT_COST_PLEVAL(hgraph, place, BBoxWChengPlEval)
    case BBoxWRSMT:
      PRINT_COST_PLEVAL(hgraph, place, BBoxWRSMTPlEval)
    case HBBox:
    case HBBoxWCheng:
    case HBBoxWRSMT:
    case HBBox0:
    case HBBox0wCheng:
    case HBBox0wRSMT:
    case HBBoxN:
    case HBBoxNwCheng:
    case HBBoxNwRSMT:
    case HBBox6:
    case HBBox6wRSMT:
      abkfatal(0, "Require regions list and assignment.");
    default:
      abkfatal(0, " Evaluator requested unknown to PlaceEvalType class\n");
  }
}

#define PRINT_COST_WREGS_PLEVAL(hgraph, regions, assignment, EVAL) \
  {                                                                \
    ::EVAL eval(hgraph, regions, assignment);                      \
    cout << eval;                                                  \
    break;                                                         \
  }

void PlaceEvalType ::costPerNet(const HGraphFixed& hgraph,
                                const LayoutBBoxes& regions,
                                vector<unsigned>& assignment) {
  switch (_type) {
    case BBoxPlEval:
    case BBoxPlEvalWt:
    case BBoxWCheng:
    case BBoxWRSMT:
      abkfatal(0, "Require placement.");
    case HBBox:
      PRINT_COST_WREGS_PLEVAL(hgraph, regions, assignment, HBBoxPlEval)
    case HBBoxWCheng:
      PRINT_COST_WREGS_PLEVAL(hgraph, regions, assignment, HBBoxWChengPlEval)
    case HBBoxWRSMT:
      PRINT_COST_WREGS_PLEVAL(hgraph, regions, assignment, HBBoxWRSMTPlEval)
    case HBBox0:
      PRINT_COST_WREGS_PLEVAL(hgraph, regions, assignment, HBBox0PlEval)
    case HBBox0wCheng:
      PRINT_COST_WREGS_PLEVAL(hgraph, regions, assignment, HBBox0wChengPlEval)
    case HBBox0wRSMT:
      PRINT_COST_WREGS_PLEVAL(hgraph, regions, assignment, HBBox0wRSMTPlEval)
    case HBBoxN:
      PRINT_COST_WREGS_PLEVAL(hgraph, regions, assignment, HBBoxNPlEval)
    case HBBoxNwCheng:
      PRINT_COST_WREGS_PLEVAL(hgraph, regions, assignment, HBBoxNwChengPlEval)
    case HBBoxNwRSMT:
      PRINT_COST_WREGS_PLEVAL(hgraph, regions, assignment, HBBoxNwRSMTPlEval)
    case HBBox6:
      PRINT_COST_WREGS_PLEVAL(hgraph, regions, assignment, HBBox6PlEval)
    case HBBox6wRSMT:
      PRINT_COST_WREGS_PLEVAL(hgraph, regions, assignment, HBBox6wRSMTPlEval)
    default:
      abkfatal(0, " Evaluator requested unknown to PlaceEvalType class\n");
  }
}

#define PRINT_COST_WRPL_PLEVAL(hgraph, regions, place, EVAL) \
  {                                                          \
    ::EVAL eval(hgraph, regions, place);                     \
    cout << eval;                                            \
    break;                                                   \
  }

void PlaceEvalType ::costPerNet(const HGraphFixed& hgraph,
                                const LayoutBBoxes& regions,
                                const Placement& place) {
  switch (_type) {
    case BBoxPlEval:
      PRINT_COST_WRPL_PLEVAL(hgraph, regions, place, BBoxPlEval)
    case BBoxPlEvalWt:
      PRINT_COST_WRPL_PLEVAL(hgraph, regions, place, BBoxPlEvalWt)
    case BBoxWCheng:
      PRINT_COST_WRPL_PLEVAL(hgraph, regions, place, BBoxWChengPlEval)
    case BBoxWRSMT:
      PRINT_COST_WRPL_PLEVAL(hgraph, regions, place, BBoxWRSMTPlEval)
    case HBBox:
      PRINT_COST_WRPL_PLEVAL(hgraph, regions, place, HBBoxPlEval)
    case HBBoxWCheng:
      PRINT_COST_WRPL_PLEVAL(hgraph, regions, place, HBBoxWChengPlEval)
    case HBBoxWRSMT:
      PRINT_COST_WRPL_PLEVAL(hgraph, regions, place, HBBoxWRSMTPlEval)
    case HBBox0:
      PRINT_COST_WRPL_PLEVAL(hgraph, regions, place, HBBox0PlEval)
    case HBBox0wCheng:
      PRINT_COST_WRPL_PLEVAL(hgraph, regions, place, HBBox0wChengPlEval)
    case HBBox0wRSMT:
      PRINT_COST_WRPL_PLEVAL(hgraph, regions, place, HBBox0wRSMTPlEval)
    case HBBoxN:
      PRINT_COST_WRPL_PLEVAL(hgraph, regions, place, HBBoxNPlEval)
    case HBBoxNwCheng:
      PRINT_COST_WRPL_PLEVAL(hgraph, regions, place, HBBoxNwChengPlEval)
    case HBBoxNwRSMT:
      PRINT_COST_WRPL_PLEVAL(hgraph, regions, place, HBBoxNwRSMTPlEval)
    case HBBox6:
      PRINT_COST_WRPL_PLEVAL(hgraph, regions, place, HBBox6PlEval)
    case HBBox6wRSMT:
      PRINT_COST_WRPL_PLEVAL(hgraph, regions, place, HBBox6wRSMTPlEval)
    default:
      abkfatal(0, " Evaluator requested unknown to PlaceEvalType class\n");
  }
}

#define COMPARISON(hgraph, regions, place, fileName, EVAL) \
  {                                                        \
    ::EVAL eval(hgraph, regions, place);                   \
    cout << eval;                                          \
    break;                                                 \
  }

void PlaceEvalType ::comparison(const HGraphFixed& hgraph,
                                const LayoutBBoxes& regions,
                                const Placement& place, const char* fileName) {
  static_cast<void>(fileName);
  switch (_type) {
    case BBoxPlEval:
      PRINT_COST_WRPL_PLEVAL(hgraph, regions, place, BBoxPlEval)
    case BBoxPlEvalWt:
      PRINT_COST_WRPL_PLEVAL(hgraph, regions, place, BBoxPlEvalWt)
    case BBoxWCheng:
      PRINT_COST_WRPL_PLEVAL(hgraph, regions, place, BBoxWChengPlEval)
    case BBoxWRSMT:
      PRINT_COST_WRPL_PLEVAL(hgraph, regions, place, BBoxWRSMTPlEval)
    case HBBox:
      PRINT_COST_WRPL_PLEVAL(hgraph, regions, place, HBBoxPlEval)
    case HBBoxWCheng:
      PRINT_COST_WRPL_PLEVAL(hgraph, regions, place, HBBoxWChengPlEval)
    case HBBoxWRSMT:
      PRINT_COST_WRPL_PLEVAL(hgraph, regions, place, HBBoxWRSMTPlEval)
    case HBBox0:
      PRINT_COST_WRPL_PLEVAL(hgraph, regions, place, HBBox0PlEval)
    case HBBox0wCheng:
      PRINT_COST_WRPL_PLEVAL(hgraph, regions, place, HBBox0wChengPlEval)
    case HBBox0wRSMT:
      PRINT_COST_WRPL_PLEVAL(hgraph, regions, place, HBBox0wRSMTPlEval)
    case HBBoxN:
      PRINT_COST_WRPL_PLEVAL(hgraph, regions, place, HBBoxNPlEval)
    case HBBoxNwCheng:
      PRINT_COST_WRPL_PLEVAL(hgraph, regions, place, HBBoxNwChengPlEval)
    case HBBoxNwRSMT:
      PRINT_COST_WRPL_PLEVAL(hgraph, regions, place, HBBoxNwRSMTPlEval)
    case HBBox6:
      PRINT_COST_WRPL_PLEVAL(hgraph, regions, place, HBBox6PlEval)
    case HBBox6wRSMT:
      PRINT_COST_WRPL_PLEVAL(hgraph, regions, place, HBBox6wRSMTPlEval)
    default:
      abkfatal(0, " Evaluator requested unknown to PlaceEvalType class\n");
  }
}

PlaceEvalXFace* PlaceEvalType::init(const HGraphFixed& hgraph,
                                    const LayoutBBoxes& regions,
                                    const Placement& place) {
  PlaceEvalXFace* plEval = NULL;
  switch (_type) {
    case BBoxPlEval:
      plEval = new ::BBoxPlEval(hgraph, regions, place);
      break;
    case BBoxPlEvalWt:
      plEval = new ::BBoxPlEvalWt(hgraph, regions, place);
      break;
    case BBoxWCheng:
      plEval = new ::BBoxWChengPlEval(hgraph, regions, place);
      break;
    case BBoxWRSMT:
      plEval = new ::BBoxWRSMTPlEval(hgraph, regions, place);
      break;
    case HBBox:
      plEval = new ::HBBoxPlEval(hgraph, regions, place);
      break;
    case HBBoxWCheng:
      plEval = new ::HBBoxWChengPlEval(hgraph, regions, place);
      break;
    case HBBoxWRSMT:
      plEval = new ::HBBoxWRSMTPlEval(hgraph, regions, place);
      break;
    case HBBox0:
      plEval = new ::HBBox0PlEval(hgraph, regions, place);
      break;
    case HBBox0wCheng:
      plEval = new ::HBBox0wChengPlEval(hgraph, regions, place);
      break;
    case HBBox0wRSMT:
      plEval = new ::HBBox0wChengPlEval(hgraph, regions, place);
      break;
    case HBBoxN:
      plEval = new ::HBBoxNPlEval(hgraph, regions, place);
      break;
    case HBBoxNwCheng:
      plEval = new ::HBBoxNwChengPlEval(hgraph, regions, place);
      break;
    case HBBoxNwRSMT:
      plEval = new ::HBBoxNwRSMTPlEval(hgraph, regions, place);
      break;
    case HBBox6:
      plEval = new ::HBBox6PlEval(hgraph, regions, place);
      break;
    case HBBox6wRSMT:
      plEval = new ::HBBox6wRSMTPlEval(hgraph, regions, place);
      break;
    default:
      abkfatal(0, " Evaluator requested unknown to PlaceEvalType class\n");
  }
  return plEval;
}

char PlaceEvalType::strbuf[32];

PlaceEvalType::operator const char*() const {
  switch (_type) {
    case BBoxPlEval:
      strncpy(strbuf, "BBoxPlEval      ", 30);
      return strbuf;
    case BBoxPlEvalWt:
      strncpy(strbuf, "BBoxPlEvalWt    ", 30);
      return strbuf;
    case BBoxWCheng:
      strncpy(strbuf, "BBoxWCheng      ", 30);
      return strbuf;
    case BBoxWRSMT:
      strncpy(strbuf, "BBoxWRSMT       ", 30);
      return strbuf;
    case HBBox:
      strncpy(strbuf, "HBBox           ", 30);
      return strbuf;
    case HBBoxWCheng:
      strncpy(strbuf, "HBBoxWCheng     ", 30);
      return strbuf;
    case HBBoxWRSMT:
      strncpy(strbuf, "HBBoxWRSMT      ", 30);
      return strbuf;
    case HBBox0:
      strncpy(strbuf, "HBBox0          ", 30);
      return strbuf;
    case HBBox0wCheng:
      strncpy(strbuf, "HBBox0wCheng    ", 30);
      return strbuf;
    case HBBox0wRSMT:
      strncpy(strbuf, "HBBox0wRSMT     ", 30);
      return strbuf;
    case HBBoxN:
      strncpy(strbuf, "HBBoxN          ", 30);
      return strbuf;
    case HBBoxNwCheng:
      strncpy(strbuf, "HBBoxNwCheng    ", 30);
      return strbuf;
    case HBBoxNwRSMT:
      strncpy(strbuf, "HBBoxNwRSMT     ", 30);
      return strbuf;
    case HBBox6:
      strncpy(strbuf, "HBBox6          ", 30);
      return strbuf;
    case HBBox6wRSMT:
      strncpy(strbuf, "HBBox6wRSMT     ", 30);
      return strbuf;
    default:
      abkfatal(0, " Evaluator requested unknown to PlaceEvalType class\n");
  }
  return strbuf;
}

ostream& operator<<(ostream& os, const PlaceEvalType& pet) {
  return os << static_cast<const char*>(pet);
}
