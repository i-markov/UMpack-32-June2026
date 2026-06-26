/**************************************************************************
***
*** Copyright (c) 2026 IonQ, Inc.
*** Copyright (c) 1995-2000 Regents of the University of California,
***               Andrew E. Caldwell, Andrew B. Kahng and Igor L. Markov
*** Copyright (c) 2000-2006 Regents of the University of Michigan,
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

// Resource-map comparison driver.
// It can build either the ISPD04 congestion map or the ISPD06 contest-style
// density map from the same placement and exposes both plotting interfaces.

#include "DB/DB.h"
#include "ISPD04CongMap.h"
#include "ISPD06DensityMap.h"
#include "ABKCommon/paramproc.h"
#include "RBPlace/RBPlacement.h"
#include "baseGeneric2DResourceMap.h"

using std::cout;
using std::endl;
using std::string;

int main(int argc, const char* argv[]) {
  NoParams noParams(argc, argv);
  BoolParam helpRequest("help", argc, argv);
  BoolParam helpRequest1("h", argc, argv);
  StringParam plotCongMap("plotCongMap", argc, argv);
  StringParam plotCongMapGP("plotCongMapGP", argc, argv);
  StringParam plotCongMapXPM("plotCongMapXPM", argc, argv);
  BoolParam noProbability("noProb", argc, argv);
  DoubleParam gridSize("gridSize", argc, argv);

  cout
      << "/********************************************************************"
         "******\n"
      << "***\n"
      << "*** Copyright (c) 1995-2000 Regents of the University of "
         "California,\n"
      << "***               Andrew E. Caldwell, Andrew B. Kahng and Igor L. "
         "Markov\n"
      << "*** Copyright (c) 2000-2006 Regents of the University of Michigan,\n"
      << "***               Saurabh N. Adya, Jarrod A. Roy, David A. Papa and\n"
      << "***               Igor L. Markov\n"
      << "***\n"
      << "***  Contact author(s): abk@cs.ucsd.edu, imarkov@umich.edu\n"
      << "***  Original Affiliation:   UCLA, Computer Science Department,\n"
      << "***                          Los Angeles, CA 90095-1596 USA\n"
      << "***\n"
      << "***  Permission is hereby granted, free of charge, to any person "
         "obtaining \n"
      << "***  a copy of this software and associated documentation files "
         "(the\n"
      << "***  Software, to deal in the Software without restriction, "
         "including\n"
      << "***  without limitation \n"
      << "***  the rights to use, copy, modify, merge, publish, distribute, "
         "sublicense, \n"
      << "***  and/or sell copies of the Software, and to permit persons to "
         "whom the \n"
      << "***  Software is furnished to do so, subject to the following "
         "conditions:\n"
      << "***\n"
      << "***  The above copyright notice and this permission notice shall be "
         "included\n"
      << "***  in all copies or substantial portions of the Software.\n"
      << "***\n"
      << "*** THE SOFTWARE IS PROVIDED AS IS, WITHOUT WARRANTY OF ANY KIND, \n"
      << "*** EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES\n"
      << "*** OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND "
         "NONINFRINGEMENT. \n"
      << "*** IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR "
         "ANY\n"
      << "*** CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF "
         "CONTRACT, TORT\n"
      << "*** OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE "
         "SOFTWARE OR\n"
      << "*** THE USE OR OTHER DEALINGS IN THE SOFTWARE.\n"
      << "***\n"
      << "***\n"
      << "*********************************************************************"
         "******\n"
      << endl;

  if (noParams.found() || helpRequest.found() || helpRequest1.found()) {
    cout << "  Congestion Map plotting utility" << endl;
    cout << "  Use '-help' or '-f filename.aux' " << endl;
    cout << "  Other options: \n"
         << "  -ispd06contest (plot overfullness as defined in the ISPD06 "
            "placement contest)\n"
         << "  -ispd04  (plot congestion as defined by Westra et al. in "
            "ISPD04)\n"
         << "  -gridXSize <unsigned> (number of grid cells in x direction)\n"
         << "  -gridYSize <unsigned> (number of grid cells in y direction)\n"
         << "  -plotMatlab base_filename  (plot congestion map in matlab "
            "format) \n"
         << "  -plotGP     base_filename (plot congestion map in gnuplot "
            "format)\n"
         << "  -plotXPM    base_filename (plot congestion map in xpm format)\n"
         << "  -targetDensity <double> (for the ISPD06 placement contest "
            "metric 0.0 - 1.0)\n"
         << "  -maxRes <double> (maximum resource usage to be used for "
            "normalization)\n"
         << endl;
    return 0;
  }
  StringParam auxFile("f", argc, argv);
  abkfatal(auxFile.found(), "must have -f <auxfilename>");

  unsigned xCells = 20;
  unsigned yCells = 20;

  UnsignedParam gridXSize("gridXSize", argc, argv);
  UnsignedParam gridYSize("gridYSize", argc, argv);

  if (gridXSize.found()) xCells = gridXSize;
  if (gridYSize.found()) yCells = gridYSize;

  BoolParam ispd04_("ispd04", argc, argv);
  BoolParam ispd06contest_("ispd06contest", argc, argv);

  bool ispd04 = false, ispd06contest = false;
  if (ispd04_.found()) ispd04 = true;
  if (ispd06contest_.found()) ispd06contest = true;
  unsigned numMaps = ispd04 ? 1 : 0 + ispd06contest ? 1 : 0;
  abkfatal(numMaps >= 1,
           "Must specify at least one resource map type with "
           "-{ispd04,ispd06contest}");
  abkfatal(numMaps >= 1,
           "Must specify at only one resource map type with "
           "-{ispd04,ispd06contest}");

  BaseGeneric2DResourceMap* resourceMap = 0;
  RBPlacement* rbplace = 0;
  DB* db = 0;

  // BoolParam check("checkDB",argc,argv);
  // DB::Parameters dbParams;
  // dbParams.alwaysCheckConsistency = check;

  // db = new DB(auxFile,dbParams);

  RBPlaceParams rbParams(argc, argv);
  rbplace = new RBPlacement(auxFile, rbParams);
  if (rbplace->getNumMacros() > 0) {
    rbplace->removeSitesBelowMacros();
    rbplace->alignCellsToRows();
  }

  if (ispd04) {
    ISPD04CongMap* ispd04map =
        new ISPD04CongMap(*rbplace, rbplace->getPlacement(), xCells, yCells);

    ispd04map->calculateDemand(rbplace->getHGraph());

    resourceMap = ispd04map;
  } else if (ispd06contest) {
    BoolParam special("specialPenalty", argc, argv);
    ISPD06DensityMap* ispd06map =
        new ISPD06DensityMap(*rbplace, special.found());

    DoubleParam targetDensity_("targetDensity", argc, argv);
    double targetDensity = 0.6;
    if (targetDensity_.found()) targetDensity = targetDensity_;
    ispd06map->setTargetDensity(targetDensity);

    resourceMap = ispd06map;
  }
  resourceMap->setNumTiles(xCells, yCells);
  resourceMap->compute();

  cout << endl << endl;
  resourceMap->printInfo(cout);
  cout << endl << endl;

  StringParam plotMatlab_("plotMatlab", argc, argv);
  if (plotMatlab_.found()) {
    resourceMap->plotMatlab(string(plotMatlab_));
  }

  StringParam plotGP_("plotGP", argc, argv);
  if (plotGP_.found()) {
    resourceMap->plotGP(string(plotGP_));
  }

  StringParam plotXPM_("plotXPM", argc, argv);
  if (plotXPM_.found()) {
    resourceMap->plotXPM(string(plotXPM_));
  }

  delete db;
  delete rbplace;
  delete resourceMap;

  return 0;
}
