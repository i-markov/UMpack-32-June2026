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

#include "baseSmallPl.h"

#include "ABKCommon/paramproc.h"

using std::cout;
using std::endl;
using std::ostream;

double BaseSmallPlacer::calculateWL(const Placement& placement) {
  double wl = 0;
  const SmallPlacementNetlist& netlist = _problem.getNetlist();

  BBox netBBox;

  for (unsigned n = netlist.getNumNets(); n != 0;) {
    n--;
    netBBox = _problem.getNetTerminalBBox(n);
    for (unsigned c = 0; c < _problem.getNumCells(); c++) {
      const BBox& pinOffset = netlist.getPinOffset(c, n);
      if (!pinOffset.isEmpty()) {
        netBBox += Point(placement[c].x + pinOffset.xMin,
                         placement[c].y + pinOffset.yMin);
        netBBox += Point(placement[c].x + pinOffset.xMax,
                         placement[c].y + pinOffset.yMax);
      }
      /**
                  if(netlist.getPinOffset(c,n).x != DBL_MAX)
                    netBBox += Point(placement[c].x+netlist.getPinOffset(c,n).x,
                                     placement[c].y+netlist.getPinOffset(c,n).y);
       **/
    }
    wl += netBBox.getHalfPerim();
  }
  return wl;
}

BaseSmallPlacer::Parameters::Parameters(int argc, const char* argv[])
    : verb(argc, argv), useInitSoln(true), minWL(None) {
  BoolParam h1("h", argc, argv);
  BoolParam h2("help", argc, argv);

  BoolParam skipInitSoln_("smPlSkipInitSoln", argc, argv);
  StringParam algo_("smPlAlgo", argc, argv);

  StringParam minWL_("smPlMinWL", argc, argv);
  BoolParam assignMethod_("smPlLeftOnly", argc, argv);
  BoolParam subsetBound_("subsetBound", argc, argv);

  StringParam startSoln_("smPlStart", argc, argv);
  StringParam splitType_("smPlSplit", argc, argv);

  if (h1.found() || h2.found()) {
    cout << "SmallPlacer Parameters" << endl;
    cout << "  -smPlSkipInitSoln" << endl;
    cout << "  -smPlAlgo {Branch | DynamicP}" << endl;

    cout << "  -smPlMinWL {None | HalfRowHeight | RowHeight | TwiceRowHeight}"
         << endl;
    cout << "  -subsetBound <bool> " << endl;

    cout << "  -smPlStart {Random | Analytic}" << endl;
    cout << "  -smPlSplit {Random | Even | All}" << endl;

    cout << "  -plot <filename> " << endl;

    cout << endl;
    return;
  }

  if (skipInitSoln_) useInitSoln = false;

  algo = BranchAndBound;
  if (algo_.found()) {
    if (strcasecmp(algo_, "Branch") == 0)
      algo = BranchAndBound;
    else if (strcasecmp(algo_, "DynamicP") == 0)
      algo = DynamicProgram;
    else
      abkfatal(0, "invalid opt algorithm specified");
  }

  assignMethod = LeftRight;
  if (assignMethod_.found()) assignMethod = LeftOnly;

  subsetBound = false;
  if (subsetBound_.found()) {
    subsetBound = true;
    abkwarn(assignMethod == LeftOnly,
            " cannot use subset bound with bidirectional assignment\n"
            " will use unidirectional assignment now\n");
    assignMethod = LeftOnly;
  }

  if (minWL_.found()) {
    if (strcasecmp(minWL_, "None") == 0)
      minWL = None;
    else if (strcasecmp(minWL_, "HalfRowHeight") == 0)
      minWL = HalfRowHeight;
    else if (strcasecmp(minWL_, "RowHeight") == 0)
      minWL = RowHeight;
    else if (strcasecmp(minWL_, "TwiceRowHeight") == 0)
      minWL = TwiceRowHeight;
    else
      abkfatal(0, "invalid spMinWL value");
  }

  splitType = AllSplit;
  if (splitType_.found()) {
    if (strcasecmp(splitType_, "Random") == 0)
      splitType = RandomSplit;
    else if (strcasecmp(splitType_, "Even") == 0)
      splitType = EvenSplit;
    else if (strcasecmp(splitType_, "All") == 0)
      splitType = AllSplit;
    else
      abkfatal(0, "invalid splitType value");
  }

  startSoln = SmartStart;
  if (startSoln_.found()) {
    if (strcasecmp(startSoln_, "Random") == 0)
      startSoln = RandomStart;
    else if (strcasecmp(startSoln_, "Analytic") == 0)
      startSoln = SmartStart;
    else
      abkfatal(0, "invalid startSoln value");
  }
}

ostream& operator<<(ostream& os, const SmallPlParams& params) {
  const char* tfStr[2] = {"false", "true"};

  os << "  useInitSoln       : " << tfStr[params.useInitSoln] << endl;

  if (params.algo == BaseSmallPlacer::Parameters::BranchAndBound) {
    os << "  assignment method : "
       << (params.assignMethod == BaseSmallPlacer::Parameters::LeftOnly
               ? "LeftOnly"
               : "LeftRight")
       << endl;
    os << "  use subset bound  : " << tfStr[params.subsetBound] << endl;
    os << "  minWirelength     : ";

    if (params.minWL == BaseSmallPlacer::Parameters::None)
      os << "None" << endl;
    else if (params.minWL == BaseSmallPlacer::Parameters::HalfRowHeight)
      os << "HalfRowHeight" << endl;
    else if (params.minWL == BaseSmallPlacer::Parameters::RowHeight)
      os << "RowHeight" << endl;
    else if (params.minWL == BaseSmallPlacer::Parameters::TwiceRowHeight)
      os << "TwiceRowHeight" << endl;
    else
      os << "Error - invalid min WL value" << endl;
  } else if (params.algo == BaseSmallPlacer::Parameters::DynamicProgram) {
    os << "  initSolnMethod    : ";
    if (params.startSoln == BaseSmallPlacer::Parameters::RandomStart)
      os << "Random" << endl;
    else if (params.startSoln == BaseSmallPlacer::Parameters::SmartStart)
      os << "Analytic" << endl;
    else
      os << "Error - invalid startSoln value" << endl;

    os << "  splitMethod       : ";
    if (params.splitType == BaseSmallPlacer::Parameters::RandomSplit)
      os << "Random" << endl;
    else if (params.splitType == BaseSmallPlacer::Parameters::EvenSplit)
      os << "Even" << endl;
    else if (params.splitType == BaseSmallPlacer::Parameters::AllSplit)
      os << "All" << endl;
    else
      os << "Error - invalid splitType value" << endl;
  } else {
    os << "Error - invalid opt algorithm value" << endl;
  }
  return os;
}
