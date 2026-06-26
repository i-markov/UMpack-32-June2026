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

// Created by Igor Markov 980320

// CHANGES
// 980325  ilm  independent fix and blk files, AUX file I/O, stream I/O
// 980402  ilm  .are files are optional now
#include <ABKCommon/abkio.h>
#include <ABKCommon/abktimer.h>
#include <ABKCommon/abktypes.h>
#include <ABKCommon/infolines.h>
#include <ABKCommon/pathDelims.h>
#include <Partitioning/partProb.h>
#include <Partitioning/partitionData.h>

#include <fstream>
#include <iomanip>
#include <cerrno>
#include <climits>
#include <string>
#include <cstdlib>

using std::cout;
using std::endl;
using std::ifstream;
using std::istream;
using std::ofstream;
using std::setw;
using std::string;
using std::vector;

static bit_vector relativeCaps;
static vector<double> tolerances;
static vector<unsigned> tolType;

static bool isAbsolutePath(const char* fileName) {
  return fileName != NULL &&
         (fileName[0] == pathDelimWindows || fileName[0] == pathDelimUnix);
}

static string extendFileName(const char* baseFileName, const char* suffix) {
  return string(baseFileName) + suffix;
}

static string makeAuxRelativeFileName(const string& dir,
                                      const char* fileName) {
  if (fileName == NULL) return string();
  if (isAbsolutePath(fileName)) return string(fileName);
  return dir + fileName;
}

static char* pathArg(string& fileName) {
  if (fileName.empty()) return NULL;
  return _abk_cpd(&fileName[0]);
}

static unsigned parseUnsignedSuffix(const char* token, unsigned prefixLen,
                                    int lineNo, const char* description) {
  abkfatal2(token != NULL && strlen(token) > prefixLen, description, lineNo);
  errno = 0;
  char* end = NULL;
  const unsigned long parsed = strtoul(token + prefixLen, &end, 10);
  abkfatal2(end != token + prefixLen && *end == '\0' && errno != ERANGE &&
                parsed <= UINT_MAX,
            description, lineNo);
  return static_cast<unsigned>(parsed);
}

PartitioningProblem ::PartitioningProblem(const char* baseFileName,
                                          Parameters params)
    : _params(params) {
  Timer tm;
  initToNULL();
  read(baseFileName);
  _postProcess();
  tm.stop();

  if (_params.verbosity.getForSysRes())
    cout << "Init+Read+PostProcess for PartProblem took " << tm << endl;
}

PartitioningProblem::PartitioningProblem(
    const char* netDFileName, const char* areFileName, const char* blkFileName,
    const char* fixFileName, const char* solFileName, Parameters params)
    : _params(params) {
  Timer tm;
  initToNULL();
  bool success = readHG(netDFileName, areFileName);
  success |= read(blkFileName, fixFileName, solFileName);

  abkwarn(success, "Reading PartitioningProblem from files may have failed");
  _postProcess();
  tm.stop();
  if (_params.verbosity.getForSysRes())
    cout << "Init+Read+PostProcess for PartProblem took " << tm << endl;
}

bool PartitioningProblem::read(const char* baseFileName) {
  Timer tm0;
  ifstream auxFile(baseFileName);
  bool success = false;
  char* newNetD = NULL;
  char* newAre = NULL;

  char* newNodes = NULL;
  char* newNets = NULL;
  char* newWts = NULL;

  char* newBlk = NULL;
  char* newFix = NULL;
  char* newSol = NULL;

  string dir;
  char* auxFileDirEnd = strrchr(const_cast<char*>(baseFileName), pathDelim);
  if (auxFileDirEnd) {
    dir.assign(baseFileName, auxFileDirEnd - baseFileName + 1);
  }

  if (auxFile) {
    bool found = false;
    int lineNo = 1;
    cout << "Reading " << baseFileName << endl;
    char oneLine[1023];
    string word1, word2;

    while (!found && !auxFile.eof()) {
      auxFile >> eathash(lineNo) >> word1 >> noeol(lineNo) >> word2 >>
          noeol(lineNo);
      abkfatal(word2 == ":",
               " Error in aux file: space-separated column expected");
      if (!strcasecmp(word1.c_str(), "CD")) {
        auxFile >> word1;
        auxFile >> needeol(lineNo++);
        if (isAbsolutePath(word1.c_str()))
          dir = word1;
        else
          dir += word1;
        if (word1[word1.length() - 1] != pathDelimWindows ||
            word1[0] == pathDelimUnix)
          dir += pathDelim;
      } else if (!strcasecmp(word1.c_str(), "PartProb")) {
        found = true;
        auxFile.getline(oneLine, 1023);
        unsigned len = strlen(oneLine), fileNum = 0;
        bool space = true;
        char* fileNames[6] = {NULL, NULL, NULL, NULL, NULL, NULL};
        unsigned j;
        for (j = 0; j != len; j++) {
          if (isspace(oneLine[j])) {
            if (!space) oneLine[j] = '\0';
            space = true;
          } else if (space) {
            space = false;
            abkfatal(fileNum < 6, " Too many filenames in AUX file");
            fileNames[fileNum++] = oneLine + j;
          }
        }

        for (j = 0; j != fileNum; j++) {
          if (strstr(fileNames[j], ".nodes")) {
            newNodes = fileNames[j];
          } else if (strstr(fileNames[j], ".NODES")) {
            newNodes = fileNames[j];
          } else if (strstr(fileNames[j], ".nets")) {
            newNets = fileNames[j];
          } else if (strstr(fileNames[j], ".NETS")) {
            newNets = fileNames[j];
          } else if (strstr(fileNames[j], ".wts")) {
            newWts = fileNames[j];
          } else if (strstr(fileNames[j], ".WTS")) {
            newWts = fileNames[j];
          }

          else if (strstr(fileNames[j], ".net")) {
            newNetD = fileNames[j];
          } else if (strstr(fileNames[j], ".NET")) {
            newNetD = fileNames[j];
          } else if (strstr(fileNames[j], ".are")) {
            newAre = fileNames[j];
          } else if (strstr(fileNames[j], ".ARE")) {
            newAre = fileNames[j];
          }

          else if (strstr(fileNames[j], ".blk")) {
            newBlk = fileNames[j];
          } else if (strstr(fileNames[j], ".BLK")) {
            newBlk = fileNames[j];
          } else if (strstr(fileNames[j], ".fix")) {
            newFix = fileNames[j];
          } else if (strstr(fileNames[j], ".FIX")) {
            newFix = fileNames[j];
          } else if (strstr(fileNames[j], ".sol")) {
            newSol = fileNames[j];
          } else if (strstr(fileNames[j], ".SOL")) {
            newSol = fileNames[j];
          }
        }

        string tmpNet = makeAuxRelativeFileName(dir, newNetD);
        string tmpAre = makeAuxRelativeFileName(dir, newAre);
        string tmpNodes = makeAuxRelativeFileName(dir, newNodes);
        string tmpNets = makeAuxRelativeFileName(dir, newNets);
        string tmpWts = makeAuxRelativeFileName(dir, newWts);
        string tmpBlk = makeAuxRelativeFileName(dir, newBlk);
        string tmpFix = makeAuxRelativeFileName(dir, newFix);
        string tmpSol = makeAuxRelativeFileName(dir, newSol);

        abkfatal(tmpNet.empty() || tmpNodes.empty(),
                 "AUX file cannot contain both net[D] and Nodes files");

        if (!tmpNet.empty()) success = readHG(pathArg(tmpNet), pathArg(tmpAre));
        if (!tmpNodes.empty())
          success =
              readHG(pathArg(tmpNodes), pathArg(tmpNets), pathArg(tmpWts));

        success |= read(pathArg(tmpBlk), pathArg(tmpFix), pathArg(tmpSol));
      }
    }
  } else {
    string newNetDStr = extendFileName(baseFileName, ".netD");
    string newAreStr = extendFileName(baseFileName, ".are");
    string newNodesStr = extendFileName(baseFileName, ".nodes");
    string newNetsStr = extendFileName(baseFileName, ".nets");
    string newWtsStr = extendFileName(baseFileName, ".wts");
    string newBlkStr = extendFileName(baseFileName, ".blk");
    string newFixStr = extendFileName(baseFileName, ".fix");
    string newSolStr = extendFileName(baseFileName, ".sol");

    ifstream netDStr(newNetDStr.c_str());
    ifstream nodesStr(newNodesStr.c_str());
    abkfatal2(netDStr || nodesStr,
              " No .aux, .nodes or .netD found: ", baseFileName);

    if (netDStr) {
      netDStr.close();
      nodesStr.close();
      success = readHG(newNetDStr.c_str(), newAreStr.c_str());
    } else {
      netDStr.close();
      nodesStr.close();
      success = readHG(newNodesStr.c_str(), newNetsStr.c_str(),
                       newWtsStr.c_str());
    }

    success |= read(newBlkStr.c_str(), newFixStr.c_str(), newSolStr.c_str());
  }
  abkwarn(success, "Reading PartitioningProblem from files may have failed\n");
  tm0.stop();

  if (_params.verbosity.getForSysRes())
    cout << "Reading PartitioningProblem took " << tm0 << endl;
  return true;
}

bool PartitioningProblem::saveAsNetDAre(const char* baseFileName) const {
  _hgraph->saveAsNetDAre(baseFileName);
  cout << " Done writing " << baseFileName << ".netD and " << baseFileName
       << ".areM" << endl;

  string newAux = extendFileName(baseFileName, ".aux");
  string newBlk = extendFileName(baseFileName, ".blk");
  string newFix = extendFileName(baseFileName, ".fix");
  string newSol = extendFileName(baseFileName, ".sol");

  const char* blkFileName = newBlk.c_str();
  const char* fixFileName = newFix.c_str();
  const char* solFileName = newSol.c_str();

  if (_terminalToBlock == NULL || _terminalToBlock->empty())
    fixFileName = NULL;
  if (_partitions == NULL) {
    blkFileName = NULL;
    solFileName = NULL;
  } else if (_bestSolnNum >= _solnBuffers->size()) {
    solFileName = NULL;
  }

  bool success = save(blkFileName, fixFileName, solFileName, true);

  abkwarn(success, "Writing PartitioningProblem to files may have failed");
  if (success) {
    ofstream auxFile(newAux.c_str());
    abkfatal(auxFile, " Can't open AUX file for writing");
    auxFile << "PartProb : " << baseFileName << ".netD " << baseFileName
            << ".areM ";
    if (blkFileName) auxFile << blkFileName << " ";
    if (fixFileName) auxFile << fixFileName << " ";
    if (solFileName) auxFile << solFileName << " ";
    auxFile << endl;
  }

  return success;
}

bool PartitioningProblem::saveAsNodesNets(const char* baseFileName) const {
  _hgraph->saveAsNodesNetsWts(baseFileName);

  cout << " Done writing " << baseFileName << ".nodes " << baseFileName
       << ".nets and " << baseFileName << ".wts" << endl;

  string newAux = extendFileName(baseFileName, ".aux");
  string newBlk = extendFileName(baseFileName, ".blk");
  string newFix = extendFileName(baseFileName, ".fix");
  string newSol = extendFileName(baseFileName, ".sol");

  const char* blkFileName = newBlk.c_str();
  const char* fixFileName = newFix.c_str();
  const char* solFileName = newSol.c_str();

  if (_terminalToBlock == NULL || _terminalToBlock->empty())
    fixFileName = NULL;
  if (_partitions == NULL) {
    blkFileName = NULL;
    solFileName = NULL;
  } else if (_bestSolnNum >= _solnBuffers->size()) {
    solFileName = NULL;
  }

  bool success = save(blkFileName, fixFileName, solFileName, false);

  abkwarn(success, "Writing PartitioningProblem to files may have failed");
  if (success) {
    ofstream auxFile(newAux.c_str());
    abkfatal(auxFile, " Can't open AUX file for writing");
    auxFile << "PartProb : " << baseFileName << ".nodes " << baseFileName
            << ".nets " << baseFileName << ".wts ";
    if (blkFileName) auxFile << blkFileName << " ";
    if (fixFileName) auxFile << fixFileName << " ";
    if (solFileName) auxFile << solFileName << " ";
    auxFile << endl;
  }

  return success;
}

bool PartitioningProblem ::save(const char* blkFileName,
                                const char* fixFileName,
                                const char* solFileName,
                                bool saveAsNetD) const {
  unsigned k, j, numWeights = (*_totalWeight).size();
  vector<double> sumCapacities(numWeights, 0);

  if (_capacities == NULL) return true;

  unsigned numPart = (*_capacities).size();

  char buf[20];

  for (k = 0; k != numPart; k++) {
    for (j = 0; j != numWeights; j++) sumCapacities[j] += (*_capacities)[k][j];
  }

  if (blkFileName && _padBlocks && _partitions) {
    unsigned numPadBlocks = (*_padBlocks).size();

    cout << " Writing " << blkFileName << " ... " << endl;
    ofstream blkFile(blkFileName);
    abkfatal(blkFile, " Can't open .blk file for writing");

    blkFile << "UCLA blk 1.0 " << endl
            << TimeStamp() << User() << Platform() << endl;

    //  Write .blk file

    blkFile << "Regular partitions : " << numPart << endl
            << "Pad partitions : " << numPadBlocks << endl
            << "Relative capacities : ";
    for (k = 0; k != (*_totalWeight).size(); k++) blkFile << "yes ";
    blkFile << endl << "Capacity tolerances : ";

    const vector<double>& firstNodeCaps = (*_capacities)[0];
    const vector<double>& firstNodeMaxCaps = (*_maxCapacities)[0];
    vector<double> capTols(firstNodeCaps.size());
    for (k = 0; k != capTols.size(); k++) {
      capTols[k] = 100.0 * fabs(firstNodeMaxCaps[k] - firstNodeCaps[k]) /
                   firstNodeCaps[k];  // sumCapacities[k] ;
      blkFile << capTols[k] << "%   ";
    }
    blkFile << endl;
    for (k = 0; k != _padBlocks->size(); k++) {
      snprintf(buf, sizeof(buf), "pb%u", k);
      blkFile << setw(7) << buf << " rect ";
      const Rectangle& rect = (*_padBlocks)[k];
      if (rect.isEmpty())
        blkFile << "        Empty ";
      else {
        blkFile << setw(9) << rect.xMin << " " << setw(9) << rect.yMin << " "
                << setw(9) << rect.xMax << " " << setw(9) << rect.yMax << " ";
      }
      //    blkFile << " : ";
      blkFile << endl;
    }

    for (k = 0; k != numPart; k++) {
      snprintf(buf, sizeof(buf), "b%u", k);
      blkFile << setw(7) << buf << " rect ";
      const Rectangle& rect = (*_partitions)[k];
      if (rect.isEmpty())
        blkFile << "        Empty ";
      else {
        blkFile << setw(9) << rect.xMin << " " << setw(9) << rect.yMin << " "
                << setw(9) << rect.xMax << " " << setw(9) << rect.yMax << " ";
      }
      blkFile << " : ";
      const vector<double>& caps = (*_capacities)[k];
      for (j = 0; j != caps.size(); j++) {
        blkFile << setw(5) << 100 * caps[j] / sumCapacities[j] << " ";
        // compute cap tols and write them here if difft from capTols
      }
      blkFile << endl;
    }
  }

  //  Now write .fix file

  if (fixFileName && _terminalToBlock && !_terminalToBlock->empty()) {
    ofstream fixFile(fixFileName);
    abkfatal(fixFile, " Can't open .fix file for writing");
    cout << " Writing " << fixFileName << " ... " << endl;
    fixFile << "UCLA fix 1.0 " << endl
            << TimeStamp() << User() << Platform() << endl;

    PartitionIds freeToMove;
    freeToMove.setToAll(numPart);

    unsigned numFixedPads = 0;
    unsigned ttBlocks = _terminalToBlock->size();
    for (k = 0; k != _hgraph->getNumTerminals(); k++) {
      if ((*_fixedConstr)[k] != freeToMove ||
          (k < ttBlocks && (*_terminalToBlock)[k] != UINT_MAX))
        numFixedPads++;
    }

    unsigned numFixedNonPads = 0;
    for (k = _hgraph->getNumTerminals(); k != _hgraph->getNumNodes(); k++) {
      PartitionIds constrainedTo = (*_fixedConstr)[k];
      if (constrainedTo != freeToMove && !constrainedTo.isEmpty())
        numFixedNonPads++;
    }

    fixFile << "Regular partitions : " << (*_partitions).size() << endl
            << "Pad partitions : " << (*_padBlocks).size() << endl;
    fixFile << "Fixed Pads : " << numFixedPads << endl;
    fixFile << "Fixed NonPads : " << numFixedNonPads << endl;

    itHGFNodeGlobal n;
    for (n = _hgraph->terminalsBegin(); n != _hgraph->terminalsEnd(); n++) {
      k = (*n)->getIndex();

      unsigned pb = UINT_MAX;
      if (k < ttBlocks) pb = (*_terminalToBlock)[k];
      PartitionIds propagatedTo = (*_fixedConstr)[k];

      if (pb != UINT_MAX || propagatedTo != freeToMove) {
        if (saveAsNetD)
          fixFile << "p" << k + 1 << " : ";
        else
          fixFile << setw(10) << _hgraph->getNodeNameByIndex((*n)->getIndex())
                  << " : ";
        if (pb != UINT_MAX) fixFile << "pb" << pb << " ";
        if (propagatedTo != freeToMove)
          for (j = 0; j != 8 * sizeof(PartitionIds); j++)
            if (propagatedTo[j]) fixFile << "b" << j << " ";
        fixFile << endl;
      }
    }

    for (; n != _hgraph->nodesEnd(); n++) {
      k = (*n)->getIndex();

      PartitionIds constrainedTo = (*_fixedConstr)[k];
      if (constrainedTo != freeToMove) {
        if (saveAsNetD)
          fixFile << "p" << k + 1 << " : ";
        else
          fixFile << setw(10) << _hgraph->getNodeNameByIndex((*n)->getIndex())
                  << " : ";
        if (constrainedTo != freeToMove && !constrainedTo.isEmpty())
          for (j = 0; j != 8 * sizeof(PartitionIds); j++)
            if (constrainedTo[j]) fixFile << "b" << j << " ";
        fixFile << endl;
      }
    }
  }

  if (solFileName && _bestSolnNum < _solnBuffers->size())
    saveBestSol(solFileName);
  return true;
}

bool PartitioningProblem::saveBestSol(const char* solFileName) const {
  unsigned bestSol = getBestSolnNum();
  ofstream solFile(solFileName);
  abkfatal(solFile, " Can't open .sol file for writing");
  cout << " Writing " << solFileName << " ... " << endl;
  solFile << "UCLA sol 1.0 " << endl
          << TimeStamp() << User() << Platform() << endl;

  unsigned numPart = (*_capacities).size();
  PartitionIds freeToMove;  //, nowhere;
  freeToMove.setToAll(numPart);

  unsigned j, k, numFixedPads = 0;

  for (k = 0; k != _hgraph->getNumTerminals(); k++) {
    if ((*_solnBuffers)[bestSol][k] != freeToMove) numFixedPads++;
  }

  unsigned numFixedNonPads = 0;
  for (k = _hgraph->getNumTerminals(); k != _hgraph->getNumNodes(); k++) {
    PartitionIds constrainedTo = (*_solnBuffers)[bestSol][k];
    if (constrainedTo != freeToMove) numFixedNonPads++;
  }

  solFile << "Regular partitions : " << (*_partitions).size() << endl
          << "Pad partitions : " << (*_padBlocks).size() << endl;
  solFile << "Fixed Pads : " << numFixedPads << endl;
  solFile << "Fixed NonPads : " << numFixedNonPads << endl;

  itHGFNodeGlobal n;
  for (n = _hgraph->terminalsBegin(); n != _hgraph->terminalsEnd(); n++) {
    k = (*n)->getIndex();

    PartitionIds propagatedTo = (*_solnBuffers)[bestSol][k];
    if (propagatedTo != freeToMove) {
      solFile << setw(10) << _hgraph->getNodeNameByIndex((*n)->getIndex())
              << " : ";
      if (propagatedTo != freeToMove)
        for (j = 0; j != 8 * sizeof(PartitionIds); j++)
          if (propagatedTo[j]) solFile << "b" << j << " ";
      solFile << endl;
    }
  }

  for (; n != _hgraph->nodesEnd(); n++) {
    k = (*n)->getIndex();

    PartitionIds constrainedTo = (*_solnBuffers)[bestSol][k];
    if (constrainedTo != freeToMove) {
      solFile << setw(10) << _hgraph->getNodeNameByIndex((*n)->getIndex())
              << " : ";
      if (constrainedTo != freeToMove)
        for (j = 0; j != 8 * sizeof(PartitionIds); j++)
          if (constrainedTo[j]) solFile << "b" << j << " ";
      solFile << endl;
    }
  }
  return true;
}

bool PartitioningProblem ::readBLK(istream& blkFile) {
  abkfatal(_hgraph != 0, " _hgraph is not yet initialized in PartProb");
  abkfatal(blkFile, " Can't open .blk file for reading");

  unsigned k, j, numWeights = _hgraph->getNumWeights();

  vector<double> minNodeWeights(numWeights, DBL_MAX);
  for (itHGFNodeGlobal v = (*_hgraph).nodesBegin(); v != (*_hgraph).nodesEnd();
       ++v) {
    for (k = 0; k != numWeights; k++) {
      double w = _hgraph->getWeight((*v)->getIndex(), k);
      if (w < minNodeWeights[k]) minNodeWeights[k] = w;
    }
  }

  int lineNo = 1;
  unsigned numPart = 0, numPadBlocks = 0;
  blkFile >> eathash(lineNo) >> needcaseword("UCLA", lineNo) >>
      needcaseword("blk", lineNo) >> skiptoeol;
  lineNo++;
  blkFile >> eathash(lineNo) >> needcaseword("Regular", lineNo) >>
      needcaseword("partitions", lineNo) >> needword(":", lineNo) >>
      my_isnumber(lineNo) >> numPart >> needeol(lineNo);

  _partitions = new vector<BBox>(numPart);
  PartitionIds freeToMove;
  freeToMove.setToAll(numPart);

  if (_fixedConstr == NULL)
    _fixedConstr = new Partitioning((*_hgraph).getNumNodes(), freeToMove);

  _capacities =
      new vector<vector<double> >(numPart, vector<double>(numWeights, -1));
  _maxCapacities =
      new vector<vector<double> >(numPart, vector<double>(numWeights, -1));
  _minCapacities =
      new vector<vector<double> >(numPart, vector<double>(numWeights, -1));
  blkFile >> eathash(lineNo) >> needcaseword("Pad", lineNo) >>
      needcaseword("partitions", lineNo) >> needword(":", lineNo) >>
      my_isnumber(lineNo) >> numPadBlocks >> needeol(lineNo);

  _padBlocks = new vector<BBox>(numPadBlocks);

  blkFile >> eathash(lineNo) >> needcaseword("Relative", lineNo) >>
      needcaseword("Capacities", lineNo) >> needword(":", lineNo);

  //    bit_vector relativeCaps(numWeights);
  relativeCaps = bit_vector(numWeights);
  for (k = 0; k != numWeights; k++) {
    string yesno;
    blkFile >> noeol(lineNo) >> isword(lineNo) >> yesno;
    if (strcasecmp(yesno.c_str(), "yes") == 0 ||
        strcasecmp(yesno.c_str(), "y") == 0)
      relativeCaps[k] = true;
    else if (strcasecmp(yesno.c_str(), "no") == 0 ||
             strcasecmp(yesno.c_str(), "n") == 0)
      relativeCaps[k] = false;
    else
      abkfatal2(0, "'yes' or 'no' expected in line ", lineNo);
  }
  blkFile >> needeol(lineNo);

  blkFile >> eathash(lineNo) >> needcaseword("Capacity", lineNo) >>
      needcaseword("tolerances", lineNo) >> needword(":", lineNo);

  //    vector<double> tolerances(numWeights);
  tolerances = vector<double>(numWeights);
  //    vector<unsigned> tolType(numWeights);
  tolType = vector<unsigned>(numWeights);
  for (k = 0; k != numWeights; k++) {
    string toler;
    blkFile >> noeol(lineNo) >> my_isnumber(lineNo) >> toler;
    string::size_type suffixPos = toler.find('b');
    if (suffixPos != string::npos) {
      string tolerValue(toler, 0, suffixPos);
      // tolerance in terms of min node width
      tolType[k] = 2;
      char* qualif = NULL;
      double fact = strtod(tolerValue.c_str(), &qualif);
      abkfatal2(fact >= 1,
                " Capacity tolerance less than min node weight in line ",
                lineNo);
      tolerances[k] = fact * minNodeWeights[k];
      abkfatal2(qualif != tolerValue.c_str(), " Expected a number in line ",
                lineNo);
    } else if ((suffixPos = toler.find('%')) != string::npos) {
      string tolerValue(toler, 0, suffixPos);
      // tolerance in relative terms
      tolType[k] = 0;
      char* qualif = NULL;
      tolerances[k] = 0.01 * strtod(tolerValue.c_str(), &qualif);
      abkfatal2(qualif != tolerValue.c_str(), " Expected a number in line ",
                lineNo);
    } else {
      char* qualif = NULL;
      tolerances[k] = strtod(toler.c_str(), &qualif);
      tolType[k] = 1;
      abkfatal2(qualif != toler.c_str(), " Expected a number in line ",
                lineNo);
      abkfatal2(static_cast<string::size_type>(qualif - toler.c_str()) ==
                    toler.size(),
                " Trailing characters near tolerance in line ", lineNo);
    }
  }
  for (k = 0; k != numPart + numPadBlocks; k++) {
    string partId, rect;
    blkFile >> eathash(lineNo) >> isword(lineNo) >> partId >> noeol(lineNo) >>
        rect >> noeol(lineNo);

    abkfatal2(strcasecmp(rect.c_str(), "rect") == 0,
              "'rect' expected in line", lineNo);
    double xMin, yMin, xMax, yMax;
    blkFile >> my_isnumber(lineNo) >> xMin >> my_isnumber(lineNo) >> yMin >>
        my_isnumber(lineNo) >> xMax >> my_isnumber(lineNo) >> yMax;
    BBox box(xMin, yMin, xMax, yMax);

    if (partId[0] == 'b') {
      unsigned partNum = parseUnsignedSuffix(
          partId.c_str(), 1, lineNo, " Bad partition Id in line ");
      abkfatal2(partNum < numPart, " Bad partition Id in line ", lineNo);
      (*_partitions)[partNum] = box;
      blkFile >> needword(":", lineNo);
      for (j = 0; j != numWeights; j++) {
        blkFile >> noeol(lineNo) >> my_isnumber(lineNo) >>
            (*_capacities)[partNum][j];
      }
    } else if (partId[0] == 'p' && partId[1] == 'b') {
      unsigned padBlockNum = parseUnsignedSuffix(
          partId.c_str(), 2, lineNo,
          " Bad pad block Id (too big a number) in line ");
      abkfatal2(padBlockNum < numPadBlocks,
                " Bad pad block Id (too big a number) in line ", lineNo);
      (*_padBlocks)[padBlockNum] = box;
    } else {
      char buf[255];
      snprintf(buf, sizeof(buf), "\n Got '%.240s' ", partId.c_str());
      abkfatal3(0, "Partition Id starting with b or pb expected in line ",
                lineNo, buf);
    }
    blkFile >> needeol(lineNo);
  }
  return true;
}

bool PartitioningProblem ::readFIX(istream& fixFile) {
  abkfatal(_hgraph != 0, " _hgraph is not yet initialized in PartProb");
  abkfatal(fixFile, " Can't open .fix file for reading");

  unsigned k, numWeights = _hgraph->getNumWeights();
  unsigned numTerminals = _hgraph->getNumTerminals();

  abkfatal(_terminalToBlock == 0, " _terminalToBlock is already init'd");
  _terminalToBlock = new vector<unsigned>(numTerminals, UINT_MAX);

  int lineNo = 1;
  unsigned numPart = 0, numPadBlocks = 0, numFixedPads = 0, numFixedNonPads = 0;

  fixFile >> eathash(lineNo);
  fixFile >> needcaseword("UCLA", lineNo) >> needcaseword("fix", lineNo);
  fixFile >> skiptoeol;
  lineNo++;

  fixFile >> eathash(lineNo);
  fixFile >> needcaseword("Regular", lineNo) >>
      needcaseword("partitions", lineNo) >> needword(":", lineNo) >>
      my_isnumber(lineNo) >> numPart >> needeol(lineNo);

  abkfatal(numPart > 1, " Need at least two regular partitions");
  abkfatal(numPart <= 32, " Can't have more than 32 regular partitions");

  if (_partitions != NULL) {
    abkfatal3((*_partitions).size() == numPart,
              "Mismatch with #partitions from blk file: line ", lineNo, "\n");
  }

  if (_capacities != NULL) {
    abkfatal2((*_capacities).size() == numPart,
              "Mismatch with #partitions from blk file: line ", lineNo);
  } else
    _capacities = new vector<vector<double> >(numPart);

  fixFile >> eathash(lineNo) >> needcaseword("Pad", lineNo) >>
      needcaseword("partitions", lineNo) >> needword(":", lineNo) >>
      my_isnumber(lineNo) >> numPadBlocks >> needeol(lineNo);

  fixFile >> eathash(lineNo) >> needcaseword("Fixed", lineNo) >>
      needcaseword("Pads", lineNo) >> needword(":", lineNo) >>
      my_isnumber(lineNo) >> numFixedPads >> needeol(lineNo);
  fixFile >> eathash(lineNo) >> needcaseword("Fixed", lineNo) >>
      needcaseword("NonPads", lineNo) >> needword(":", lineNo) >>
      my_isnumber(lineNo) >> numFixedNonPads >> needeol(lineNo);

  PartitionIds freeToMove;
  freeToMove.setToAll(numPart);
  if (_fixedConstr == NULL) {
    _fixedConstr = new Partitioning((*_hgraph).getNumNodes(), freeToMove);
  }

  for (k = 0; k != numFixedPads + numFixedNonPads; k++) {
    string nodeId;
    fixFile >> eathash(lineNo) >> isword(lineNo) >> nodeId >> noeol(lineNo) >>
        needword(":", lineNo) >> noeol(lineNo);

    const HGFNode& node = _hgraph->getNodeByName(nodeId.c_str());

    if (!_hgraph->isTerminal(node.getIndex())) {
      PartitionIds constr;

      while (fixFile.peek() != '\n' && fixFile.peek() != '\r') {
        string partId;
        fixFile >> isword(lineNo) >> partId >> eatblank;
        abkfatal3(partId[0] == 'b',
                  " NonPad node can be constrained to regular partitions only ",
                  "in line ", lineNo);
        unsigned partNum = parseUnsignedSuffix(
            partId.c_str(), 1, lineNo, "Bad partition Id in line ");
        abkfatal2(partNum < numPart, "Bad partition Id in line ", lineNo);
        constr.addToPart(partNum);
      }
      if (constr.isEmpty()) constr = freeToMove;
      (*_fixedConstr)[node.getIndex()] = constr;
    } else {
      PartitionIds constr;

      while (fixFile.peek() != '\n' && fixFile.peek() != '\r') {
        string partId;
        fixFile >> isword(lineNo) >> partId >> eatblank;
        if (partId[0] == 'b') {
          unsigned partNum = parseUnsignedSuffix(
              partId.c_str(), 1, lineNo,
              "Partition Id out of range in line ");
          abkfatal2(partNum < numPart, "Partition Id out of range in line ",
                    lineNo);
          constr.addToPart(partNum);
        } else if (partId[0] == 'p' && partId[1] == 'b') {
          unsigned blockNum = parseUnsignedSuffix(
              partId.c_str(), 2, lineNo, "Block Id out of range in line ");
          abkfatal2(blockNum < numPadBlocks,
                    "Block Id out of range in line ", lineNo);
          if ((*_terminalToBlock)[node.getIndex()] != UINT_MAX) {
            char errMess[256];
            snprintf(errMess, sizeof(errMess),
                     " Pad %.180s is constrained to >1 pad block in line %d\n",
                     nodeId.c_str(), lineNo);
            abkfatal(0, errMess);
          }
          (*_terminalToBlock)[node.getIndex()] = blockNum;
          for (unsigned kw = 0; kw != numWeights; kw++)
            _hgraph->setWeight(node.getIndex(), 0.0, kw);
        } else {
          char buf[255];
          snprintf(buf, sizeof(buf), ".\n Got %.240s ", partId.c_str());
          abkfatal3(0, " Block or partition Id expected in line ", lineNo, buf);
        }
      }
      if (constr.isEmpty()) constr = freeToMove;
      (*_fixedConstr)[node.getIndex()] = constr;
    }

    fixFile >> needeol(lineNo);
  }
  return true;
}

bool PartitioningProblem ::readSOL(istream& solFile) {
  abkfatal(_hgraph != 0, " _hgraph is not yet initialized in PartProb");
  abkfatal(solFile, " Can't open .sol file for reading");

  unsigned k;

  int lineNo = 1;
  unsigned numPart = 0, numPadBlocks = 0, numFixedPads = 0, numFixedNonPads = 0;

  solFile >> eathash(lineNo);
  solFile >> needcaseword("UCLA", lineNo) >> needcaseword("sol", lineNo);
  solFile >> skiptoeol;
  lineNo++;

  solFile >> eathash(lineNo);
  solFile >> needcaseword("Regular", lineNo) >>
      needcaseword("partitions", lineNo) >> needword(":", lineNo) >>
      my_isnumber(lineNo) >> numPart >> needeol(lineNo);

  abkfatal(numPart > 1, " Need at least two regular partitions");
  abkfatal(numPart <= 32, " Can't have more than 32 regular partitions");
  PartitionIds freeToMove, nowhere;
  freeToMove.setToAll(numPart);

  if (!_solnBuffers) {
    _solnBuffers = new PartitioningBuffer(_hgraph->getNumNodes(), 1);
    unsigned kk = 0;
    for (; kk != _hgraph->getNumTerminals(); ++kk)
      (*_solnBuffers)[0][kk] = nowhere;
    for (; kk != _hgraph->getNumNodes(); ++kk)
      (*_solnBuffers)[0][kk] = freeToMove;
  }

  if (_partitions != NULL) {
    abkfatal3((*_partitions).size() == numPart,
              "Mismatch with #partitions from blk file: line ", lineNo, "\n");
  }

  if (_fixedConstr) {
    if (_fixedConstr == NULL)
      _fixedConstr = new Partitioning((*_hgraph).getNumNodes(), freeToMove);
  }

  if (_capacities != NULL) {
    abkfatal2((*_capacities).size() == numPart,
              "Mismatch with #partitions from blk file: line ", lineNo);
  } else
    _capacities = new vector<vector<double> >(numPart);

  solFile >> eathash(lineNo) >> needcaseword("Pad", lineNo) >>
      needcaseword("partitions", lineNo) >> needword(":", lineNo) >>
      my_isnumber(lineNo) >> numPadBlocks >> needeol(lineNo);

  solFile >> eathash(lineNo) >> needcaseword("Fixed", lineNo) >>
      needcaseword("Pads", lineNo) >> needword(":", lineNo) >>
      my_isnumber(lineNo) >> numFixedPads >> needeol(lineNo);
  solFile >> eathash(lineNo) >> needcaseword("Fixed", lineNo) >>
      needcaseword("NonPads", lineNo) >> needword(":", lineNo) >>
      my_isnumber(lineNo) >> numFixedNonPads >> needeol(lineNo);

  for (k = 0; k != numFixedPads + numFixedNonPads; k++) {
    string nodeId;
    solFile >> eathash(lineNo) >> isword(lineNo) >> nodeId >> noeol(lineNo) >>
        needword(":", lineNo) >> noeol(lineNo);

    PartitionIds constr;
    const HGFNode& node = _hgraph->getNodeByName(nodeId.c_str());
    if (!_hgraph->isTerminal(node.getIndex())) {
      while (solFile.peek() != '\n' && solFile.peek() != '\r') {
        string partId;
        solFile >> isword(lineNo) >> partId >> eatblank;
        abkfatal3(partId[0] == 'b',
                  " NonPad node can be constrained to regular partitions only ",
                  "in line ", lineNo);
        unsigned partNum = parseUnsignedSuffix(
            partId.c_str(), 1, lineNo, "Bad partition Id in line ");
        abkfatal2(partNum < numPart, "Bad partition Id in line ", lineNo);
        constr.addToPart(partNum);
      }
    } else {
      while (solFile.peek() != '\n' && solFile.peek() != '\r') {
        string partId;
        solFile >> isword(lineNo) >> partId >> eatblank;
        if (partId[0] == 'b') {
          unsigned partNum = parseUnsignedSuffix(
              partId.c_str(), 1, lineNo,
              "Partition Id out of range in line ");
          abkfatal2(partNum < numPart, "Partition Id out of range in line ",
                    lineNo);
          constr.addToPart(partNum);
        } else if (partId[0] == 'p' && partId[1] == 'b') {
          // ignore
        } else {
          char buf[255];
          snprintf(buf, sizeof(buf), ".\n Got %.240s ", partId.c_str());
          abkfatal3(0, " Block or partition Id expected in line ", lineNo, buf);
        }
      }
    }

    if (!constr.isEmpty())
      (*_solnBuffers)[0][node.getIndex()] = constr;
    else
      (*_solnBuffers)[0][node.getIndex()] = freeToMove;

    solFile >> needeol(lineNo);
  }
  for (k = 1; k < _solnBuffers->size(); k++)
    (*_solnBuffers)[k] = (*_solnBuffers)[0];
  return true;
}

bool PartitioningProblem::readHG(const char* f1, const char* f2,
                                 const char* f3) {
  _hgraph = new HGraphFixed(f1, f2, f3, _params);
  cout << " Done reading HGraph" << endl;
  return true;
}

bool PartitioningProblem::read(const char* blkFileName, const char* fixFileName,
                               const char* solFileName) {
  abkfatal(_hgraph != NULL, "Must call readHGraph first");
  abkfatal(blkFileName, "Missing .blk file");

  _ownsData = true;
  _costs = vector<double>(1, DBL_MAX);
  _violations = vector<double>(1, DBL_MAX);
  _imbalances = vector<double>(1, DBL_MAX);

  unsigned numWeights = _hgraph->getNumWeights();

  itHGFNodeGlobalMutable v = _hgraph->nodesBegin();

  if (blkFileName) {
    cout << " Reading " << blkFileName << " ... " << endl;
    ifstream ifs(blkFileName);
    readBLK(ifs);
  }

  if (fixFileName) {
    cout << " Reading " << fixFileName << " ... " << endl;
    ifstream ifs(fixFileName);
    readFIX(ifs);
  } else  // still need _terminalToBlock
  {
    abkwarn(_hgraph->getNumTerminals() == 0,
            "Terminals specified in .net or .netD file, but "
            "no .fix file given.\n  Making all nodes nonterminal");

    _hgraph->clearTerminals();  // set num terminals to 0
    abkfatal(_terminalToBlock == 0, " _terminalToBlock is already init'd");
    _terminalToBlock =
        new vector<unsigned>(_hgraph->getNumTerminals(), UINT_MAX);
  }
  // count the total weight of the nodes in the problem
  _totalWeight = new vector<double>(numWeights, 0.0);
  for (v = _hgraph->nodesBegin(); v != _hgraph->nodesEnd(); v++) {
    for (unsigned k = 0; k != numWeights; k++)
      (*_totalWeight)[k] += _hgraph->getWeight((*v)->getIndex(), k);
  }

  unsigned k, j, numPart = _capacities->size();

  for (k = 0; k != numWeights; k++) {
    for (j = 0; j != numPart; j++) {
      if (relativeCaps[k]) (*_capacities)[j][k] *= 0.01 * (*_totalWeight)[k];
      switch (tolType[k]) {
        case 0:
          (*_maxCapacities)[j][k] = (*_capacities)[j][k] * (1 + tolerances[k]);
          (*_minCapacities)[j][k] = (*_capacities)[j][k] * (1 - tolerances[k]);
          break;
        case 1:
        case 2:
          (*_maxCapacities)[j][k] = (*_capacities)[j][k] + tolerances[k];
          (*_minCapacities)[j][k] = (*_capacities)[j][k] - tolerances[k];
          break;
        default:
          abkfatal(0, "Unknown balance tolerance type");
      }
    }
  }
  PartitionIds freeToMove, nowhere;
  freeToMove.setToAll(numPart);
  if (_capacities) {
    _solnBuffers = new PartitioningBuffer(_hgraph->getNumNodes(), 1);
    for (k = 0; k != _hgraph->getNumTerminals(); k++)
      (*_solnBuffers)[0][k] = nowhere;
    for (; k != _hgraph->getNumNodes(); k++) (*_solnBuffers)[0][k] = freeToMove;
  }

  if (solFileName) {
    cout << " Reading " << solFileName << " ... " << endl;
    ifstream ifs(solFileName);
    readSOL(ifs);
  }

  if (!_solnBuffers) {
    _solnBuffers = new PartitioningBuffer(_hgraph->getNumNodes(), 1);
    for (k = 0; k != _hgraph->getNumTerminals(); k++)
      (*_solnBuffers)[0][k] = nowhere;
    for (; k != _hgraph->getNumNodes(); k++) (*_solnBuffers)[0][k] = freeToMove;
  }

  return true;
}
