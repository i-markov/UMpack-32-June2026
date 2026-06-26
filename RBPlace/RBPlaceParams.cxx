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

// Created 051103 by David Papa to move all parameter related code and
//       definitions out of RBPlacement.h

#include <cerrno>
#include <climits>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <set>
#include <string>
#include <strings.h>

#include "ABKCommon/paramproc.h"
#include "RBPlacement.h"

using std::endl;
using std::ostream;
using std::string;
using std::vector;

namespace {

bool parsePositiveUnsigned(const char* value, unsigned& parsed) {
  if (value == NULL || *value == '\0') return false;
  errno = 0;
  char* end = NULL;
  unsigned long result = strtoul(value, &end, 10);
  if (end == value || *end != '\0' || result == 0 || errno == ERANGE ||
      result > UINT_MAX)
    return false;
  parsed = static_cast<unsigned>(result);
  return true;
}

}  // namespace

#ifdef _MSC_VER
#pragma warning(disable : 4355)
#endif

RBPlacement::Parameters::Parameters()
    : verb("1 1 1"),
      numRowsToRemove(0),
      spaceCellsAlg(EQUAL_SPACE),
      remCongestion(0),
      adaptivePlotting(true),
      compressedPlotting(true),
      xpixels(1280),
      ypixels(1024),
      overlapCells(100000),
      greedyMoveCells(15000),
      macroOverlapDecreasing(true),
#ifdef USEFLOORIST
      enableFloorist(false),
#endif
      centerOfGravityFSFPEval(false),
      pinToPinFSFPEval(false),
      obstacleFile(""),
      allowSoftBlocks(true) {
}

RBPlacement::Parameters::Parameters(int argc, const char* argv[])
    : HGraphParameters(argc, argv),
      verb(argc, argv),
      numRowsToRemove(0),
      spaceCellsAlg(EQUAL_SPACE),
      remCongestion(0),
      adaptivePlotting(true),
      compressedPlotting(true),
      xpixels(1280),
      ypixels(1024),
      overlapCells(100000),
      greedyMoveCells(15000),
      macroOverlapDecreasing(true),
#ifdef USEFLOORIST
      enableFloorist(false),
#endif
      centerOfGravityFSFPEval(false),
      pinToPinFSFPEval(false),
      obstacleFile(""),
      allowSoftBlocks(true) {
  UnsignedParam numRowsToRemove_("numRowsToRemove", argc, argv);
  StringParam spaceCellsAlg_("spaceCellsAlg", argc, argv);
  BoolParam remCongestion_("remCongestion", argc, argv);
  BoolParam noAdaptivePlotting("noAdaptivePlotting", argc, argv);
  BoolParam noCompressedPlotting("noCompressedPlotting", argc, argv);
  BoolParam centerOfGravityFSFP("centerOfGravityFSFP", argc, argv);
  BoolParam pinToPinFSFP("pinToPinFSFP", argc, argv);

  if (remCongestion_.found()) remCongestion = 1;

  if (numRowsToRemove_.found()) numRowsToRemove = numRowsToRemove_;

  if (spaceCellsAlg_.found()) {
    if (strcasecmp(spaceCellsAlg_, "EQUAL") == 0)
      spaceCellsAlg = EQUAL_SPACE;
    else if (strcasecmp(spaceCellsAlg_, "PIN1") == 0)
      spaceCellsAlg = WITH_PIN_ALG1;
    else if (strcasecmp(spaceCellsAlg_, "PIN2") == 0)
      spaceCellsAlg = WITH_PIN_ALG2;
    else
      abkwarn(0, "Unrecognize spaceCellsAlg! Using default value");
  }

  if (noAdaptivePlotting.found()) {
    adaptivePlotting = false;
    compressedPlotting = false;
  } else {
    StringParam adaptivePlotting_("adaptivePlotting", argc, argv);
    if (adaptivePlotting_.found()) {
      const string resolution = static_cast<const char*>(adaptivePlotting_);
      const string::size_type separator =
          resolution.find_first_not_of("0123456789");
      unsigned xpixeltmp = 0;
      unsigned ypixeltmp = 0;
      if (separator == string::npos || separator + 1 == resolution.size() ||
          resolution.find_first_not_of("0123456789", separator + 1) !=
              string::npos ||
          !parsePositiveUnsigned(resolution.substr(0, separator).c_str(),
                                 xpixeltmp) ||
          !parsePositiveUnsigned(resolution.substr(separator + 1).c_str(),
                                 ypixeltmp)) {
        abkwarn(0, "Improper resolution given");
      } else {
        xpixels = xpixeltmp;
        ypixels = ypixeltmp;
      }
    }
  }

  if (noCompressedPlotting.found()) {
    compressedPlotting = false;
  }

  if (centerOfGravityFSFP.found()) centerOfGravityFSFPEval = true;
  if (pinToPinFSFP.found()) pinToPinFSFPEval = true;

  UnsignedParam overlapCells_("overlapRemCutoff", argc, argv);
  if (overlapCells_.found()) {
    overlapCells = overlapCells_;
  }

  BoolParam faster_("faster", argc, argv);
  BoolParam ispd06("ispd06", argc, argv);
  if (faster_.found() || ispd06.found()) {
    greedyMoveCells = 5000;
  }

  UnsignedParam greedyMoveCells_("greedyMoveCells", argc, argv);
  if (greedyMoveCells_.found()) {
    greedyMoveCells = greedyMoveCells_;
  }

  StringParam obstacleFile_("obstacleFile", argc, argv);
  if (obstacleFile_.found()) {
    obstacleFile = obstacleFile_;
  }

  BoolParam noSoftBlocks("noSoftBlocks", argc, argv);
  if (noSoftBlocks.found()) {
    allowSoftBlocks = false;
  }

  BoolParam macroOverlapIncreasing("macroOverlapIncreasing", argc, argv);
  if (macroOverlapIncreasing.found()) {
    macroOverlapDecreasing = false;
  }

#ifdef USEFLOORIST
  BoolParam floorist("Floorist", argc, argv);
  if (floorist.found()) {
    enableFloorist = true;
  }
#endif
}

RBPlacement::Parameters::Parameters(const Parameters& orig)
    : HGraphParameters(orig),
      verb(orig.verb),
      numRowsToRemove(orig.numRowsToRemove),
      spaceCellsAlg(orig.spaceCellsAlg),
      remCongestion(orig.remCongestion),
      adaptivePlotting(orig.adaptivePlotting),
      compressedPlotting(orig.compressedPlotting),
      xpixels(orig.xpixels),
      ypixels(orig.ypixels),
      overlapCells(orig.overlapCells),
      greedyMoveCells(orig.greedyMoveCells),
      macroOverlapDecreasing(orig.macroOverlapDecreasing),
#ifdef USEFLOORIST
      enableFloorist(orig.enableFloorist),
#endif
      centerOfGravityFSFPEval(orig.centerOfGravityFSFPEval),
      pinToPinFSFPEval(orig.pinToPinFSFPEval),
      obstacleFile(orig.obstacleFile),
      allowSoftBlocks(orig.allowSoftBlocks) {
}

ostream& operator<<(ostream& out, const RBPlaceParams& params) {
  out << " RBPlacement Parameters" << endl;
  out << "  Verbosity:         " << params.verb << endl;
  out << "  NumRowsToRemove:   " << params.numRowsToRemove << endl;
  out << "  SpaceCellsAlg:     ";
  switch (params.spaceCellsAlg) {
    case RBPlacement::Parameters::EQUAL_SPACE:
      out << "Equal Spacing\n";
      break;
    case RBPlacement::Parameters::WITH_PIN_ALG1:
      out << "With Pin 1st algorithm\n";
      break;
    case RBPlacement::Parameters::WITH_PIN_ALG2:
      out << "With Pin 2nd algorithm\n";
      break;
    default:
      out << "Unknown Spacing Alg\n";
      break;
  }
  out << "  Remove Congestion: " << params.remCongestion << endl;
  out << static_cast<HGraphParameters>(params) << endl;
  out << "  adaptivePlotting: " << (params.adaptivePlotting ? "true" : "false")
      << endl;
  if (params.adaptivePlotting) {
    out << "    X pixels: " << params.xpixels << endl;
    out << "    Y pixels: " << params.ypixels << endl;
  }
  out << "  compressedPlotting: "
      << (params.compressedPlotting ? "true" : "false") << endl;
  if (params.centerOfGravityFSFPEval || params.pinToPinFSFPEval)
    out << "  FSFP HPWL Evaluation: "
        << (params.centerOfGravityFSFPEval ? "Center of Gravity" : "Pin to Pin")
        << endl;
  out << "  Obstacle File: " << params.obstacleFile << endl;
  return out;
}
