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

// CHANGES
// ilm 293000  in addition to '\n', the parser needs to check for '\r'

#include <ABKCommon/abkassert.h>
#include <ABKCommon/abkio.h>
#include <ABKCommon/pathDelims.h>
#include <Ctainers/bitBoard.h>
#include <HGraph/hgFixed.h>

#include <cerrno>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <strings.h>

using std::cerr;
using std::cout;
using std::endl;
using std::ifstream;
using std::pair;
using std::stringstream;
using std::string;
using std::vector;

vector<pair<HGFNode*, char> > _USELESS_declaration_for_compilerBut_;
pair<HGFNode*, char> _USELESS_declaration_for_compilerBut_NUmber2_;

namespace {

std::string resolvePath(const std::string& dir, const char* fileName) {
  if (!fileName) return std::string();
  if (fileName[0] == pathDelimWindows || fileName[0] == pathDelimUnix)
    return fileName;
  return dir + fileName;
}

// Parse and validate the node index encoded in a .netD/.areM token such as
// "p123" (pad, 1-based) or "a45" (cell, 0-based offset past the terminals).
// The historical code used atoi() with no error or range checking and relied on
// getNodeByIdx()'s abkassert3 bound, which compiles to a no-op in optimized
// builds -- so a corrupt or crafted file could drive an out-of-bounds access.
// This validates the suffix (strtol + errno + end-of-token) and enforces the
// range with an always-on abkfatal.
unsigned parseNodeIndexToken(const string& token, unsigned numTerminals,
                             unsigned numNodes, const char* fileName) {
  const char* digits = token.c_str() + 1;
  errno = 0;
  char* end = 0;
  long parsed = strtol(digits, &end, 10);
  abkfatal3(end != digits && *end == '\0' && errno == 0,
            "malformed node index token in ", fileName, token);
  abkfatal3(parsed >= 0, "negative node index in ", fileName, token);

  unsigned nodeIdx;
  if (token[0] == 'p')
    nodeIdx = static_cast<unsigned>(parsed) - 1;
  else
    nodeIdx = static_cast<unsigned>(parsed) + numTerminals;

  abkfatal3(nodeIdx < numNodes, "out-of-range node index in ", fileName, token);
  return nodeIdx;
}

}  // namespace

HGraphFixed::HGraphFixed(const char* Filename1, const char* Filename2,
                         const char* Filename3, HGraphParameters param)
    : HGraphBase(param) {
  abkfatal(Filename1, "null char* for first filename in HGraph ctor");

  const char* auxFileName = 0;
  const char* netDFileName = 0;
  const char* areMFileName = 0;
  const char* wtsFileName = 0;
  const char* netsFileName = 0;
  const char* nodesFileName = 0;

  const char* files[3] = {Filename1, Filename2, Filename3};

  for (unsigned f = 0; f < 3; f++) {
    if (!files[f]) continue;

    if (strstr(files[f], ".aux"))
      auxFileName = files[f];
    else if (strstr(files[f], ".wts"))
      wtsFileName = files[f];
    else if (strstr(files[f], ".are"))
      areMFileName = files[f];
    else if (strstr(files[f], ".nodes"))
      nodesFileName = files[f];
    else if (strstr(files[f], ".netD"))
      netDFileName = files[f];
    else if (strstr(files[f], ".nets"))
      netsFileName = files[f];
    else if (strstr(files[f], ".net"))
      netDFileName = files[f];
  }

  abkfatal(!((netDFileName || areMFileName) &&
             (netsFileName || nodesFileName || wtsFileName)),
           "can't mix netD/areM and nets/nodes/wts");

  if (auxFileName)
    parseAux(auxFileName);
  else if (netDFileName) {
    readNetD(netDFileName);
    if (areMFileName)
      readAreM(areMFileName);
    else {  // set unit weights
      loadUnitWeights();
      cout << "are(M) file not found. loading unit weights for all nodes"
           << endl;
    }
  } else if (nodesFileName) {
    abkfatal(netsFileName, "must have both a nodes and nets file");

    readNodes(nodesFileName);
    readNets(netsFileName);
    if (wtsFileName)
      readWts(wtsFileName);
    else {  // set unit weights
      loadUnitWeights();
    }
  } else {
    cerr << "aux, netD or nets file not found.  ";
    cerr << "Filename arguments were:" << endl;
    if (Filename1) cout << Filename1 << endl;
    if (Filename2) cout << Filename2 << endl;
    if (Filename3) cout << Filename3 << endl;
  }
  finalize();
}

void HGraphFixed::parseAux(const char* baseFileName) {
  ifstream auxFile(baseFileName);
  char* newNetD = NULL;
  char* newAreM = NULL;
  char* newNets = NULL;
  char* newNodes = NULL;
  char* newWts = NULL;

  std::string dir;

  const char* auxFileDirEnd = strrchr(baseFileName, pathDelim);
  if (auxFileDirEnd) {
    string tmp = baseFileName;
    tmp.resize(auxFileDirEnd - baseFileName + 1);
    dir = tmp;
  }

  if (auxFile) {
    bool found = false;
    int lineno = 1;
    cout << "Reading " << baseFileName << endl;
    char oneLine[1023];
    string word1, word2;

    while (!found && !auxFile.eof()) {
      auxFile >> eathash(lineno) >> word1 >> noeol(lineno) >> word2 >>
          eatblank >> noeol(lineno);
      abkfatal2(
          word2 == ":",
          " Error in aux file: space-separated column expected. got:", word2);
      if (!strcasecmp(word1.c_str(), "CD")) {
        auxFile >> word1;
        auxFile >> needeol(lineno++);
        if (word1[0] == pathDelimWindows || word1[0] == pathDelimUnix)
          dir = word1;
        else
          dir += word1;
        char fDel[2] = {pathDelim, '\0'};
        if (word1[word1.length() - 1] != pathDelimWindows ||
            word1[0] == pathDelimUnix)
          dir += fDel;
      } else if (!strcasecmp(word1.c_str(), "HGraph") ||
                 !strcasecmp(word1.c_str(), "HGraphWPins") ||
                 !strcasecmp(word1.c_str(), "PartProb") ||
                 !strcasecmp(word1.c_str(), "FPPROBLEM") ||
                 !strcasecmp(word1.c_str(),
                             "RowBasedPlacement"))  // used in RBPl ctor
      {
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
          if (strstr(fileNames[j], ".nets"))
            newNets = fileNames[j];
          else if (strstr(fileNames[j], ".NETS"))
            newNets = fileNames[j];
          else if (strstr(fileNames[j], ".nodes"))
            newNodes = fileNames[j];
          else if (strstr(fileNames[j], ".NODES"))
            newNodes = fileNames[j];
          else if (strstr(fileNames[j], ".wts"))
            newWts = fileNames[j];
          else if (strstr(fileNames[j], ".WTS"))
            newWts = fileNames[j];
          else if (strstr(fileNames[j], ".net"))
            newNetD = fileNames[j];
          else if (strstr(fileNames[j], ".NET"))
            newNetD = fileNames[j];
          else if (strstr(fileNames[j], ".are"))
            newAreM = fileNames[j];
          else if (strstr(fileNames[j], ".ARE"))
            newAreM = fileNames[j];
        }

        abkfatal(!((newNetD || newAreM) && (newNets || newNodes || newWts)),
                 "can't mix netD/areM and nets/nodes/wts");

        std::string tmpNetD = resolvePath(dir, newNetD);
        std::string tmpAreM = resolvePath(dir, newAreM);
        std::string tmpNets = resolvePath(dir, newNets);
        std::string tmpNodes = resolvePath(dir, newNodes);
        std::string tmpWts = resolvePath(dir, newWts);

        if (!tmpNetD.empty()) {
          readNetD(tmpNetD.c_str());
          if (!tmpAreM.empty())
            readAreM(tmpAreM.c_str());
          else {  // set unit weights
            loadUnitWeights();
            cout << "are(M) file not found. loading unit "
                 << "weights for all nodes" << endl;
          }
        } else if (!tmpNodes.empty()) {
          abkassert(!tmpNets.empty(), "must have .nets file with .nodes");
          readNodes(tmpNodes.c_str());
          readNets(tmpNets.c_str());
          if (!tmpWts.empty())
            readWts(tmpWts.c_str());
          else {  // set unit weights
            loadUnitWeights();
          }
        } else
          abkfatal(0, "found neither net[D] or nodes files");
      }
    }
    cout << "Finished reading " << baseFileName << " and referenced files"
         << endl;
  } else
    abkfatal2(0, baseFileName, "  not found");
}

void HGraphFixed::readNetD(const char* netDFileName) {
  ifstream netD(netDFileName);
  abkfatal2(netD, " Could not open ", netDFileName);

  char t;
  netD >> t;  // eat 1st 0
  unsigned numNets, numModules, padOffset, numWeights = 0;
  netD >> _numPins >> numNets >> numModules >> padOffset;
  abkfatal2(netD, " malformed header in ", netDFileName);
  abkfatal2(numModules > 0, " numModules must be > 0 in ", netDFileName);
  if (padOffset == 0) padOffset = numModules - 1;
  abkfatal2(padOffset < numModules, " padOffset out of range in ", netDFileName);
  _numTerminals = numModules - padOffset - 1;

  init(numModules, numWeights);

  // Build each node name once into a reusable buffer. The previous code
  // constructed a fresh stringstream per node and called .str().c_str() twice
  // (two temporary std::strings per node) -- the dominant cost when loading a
  // multi-million-node .netD. Reserve both containers up front (matching the
  // .nodes reader).
  _nodeNames.reserve(numModules);
  _nodeNamesMap.reserve(numModules);
  char nameBuf[24];
  for (unsigned n = 0; n < numModules; n++) {
    if (n < _numTerminals)
      snprintf(nameBuf, sizeof(nameBuf), "p%u", n + 1);
    else
      snprintf(nameBuf, sizeof(nameBuf), "a%u", n - _numTerminals);

    _nodeNames.push_back(nameBuf);
    _nodeNamesMap[nameBuf] = n;
  }

  HGFEdge* edge = 0;
  int lineno = 6;
  string tmp;
  vector<pair<HGFNode*, char> > netPins;

  unsigned numEdgesSeen = 0;

  while (netD >> tmp) {
    if (tmp.empty()) break;
    netD >> t;

    t = toupper(t);
    if (t == 'S') {
      if (++numEdgesSeen > numNets) {
        abkwarn3(0, "line ", lineno, ": extra edges in .netD file");
        break;
      }

      if (netPins.size() > 1)  // at the start of a new net,
                               // construct the one that just finished
      {
        edge = addEdge();
        for (unsigned i = 0; i < netPins.size(); i++) {
          const char pinDir =
              _param.makeAllSrcSnk ? 'B' : netPins[i].second;
          switch (pinDir) {
            case 'O':
#ifdef SIGNAL_DIRECTIONS
              addSrc(*netPins[i].first, *edge);
              break;
#endif
            case 'I':
#ifdef SIGNAL_DIRECTIONS
              addSnk(*netPins[i].first, *edge);
              break;
#endif
            case 'B':
            case '1':
            case 'L':  // !!! why it happens???
              addSrcSnk(*netPins[i].first, *edge);
              break;
            default:  // error
              abkfatal3(0, t, ": unexpected symbol encountered in line ",
                        lineno);
          }
        }
      }

      netPins.clear();
    }
    netD >> t;

    char origchar = t;
    t = toupper(t);

    if (netD)  // it may be that we just finished the file
    {
      if (t != 'B' && t != 'I' && t != 'O' && t != '1') {
        netD.putback(origchar);
        abkfatal(netD, "could not putback to stream");
        t = '1';
      }
    }

    unsigned nodeIdx =
        parseNodeIndexToken(tmp, _numTerminals, getNumNodes(), netDFileName);

    HGFNode& node = getNodeByIdx(nodeIdx);

    netPins.push_back(pair<HGFNode*, char>(&node, t));
    lineno++;
  }

  if (netPins.size() > 1) {
    edge = addEdge();
    for (unsigned i = 0; i < netPins.size(); i++) {
      const char pinDir = _param.makeAllSrcSnk ? 'B' : netPins[i].second;
      switch (pinDir) {
        case 'O':
#ifdef SIGNAL_DIRECTIONS
          addSrc(*netPins[i].first, *edge);
          break;
#endif
        case 'I':
#ifdef SIGNAL_DIRECTIONS
          addSnk(*netPins[i].first, *edge);
          break;
#endif
        case 'B':
        case '1':
        case 'L':  // !!! why it happens???
          addSrcSnk(*netPins[i].first, *edge);
          break;
        default:  // error
          abkfatal3(0, t, ": unexpected symbol encountered in line ", lineno);
      }
    }
  }
  _netNames = vector<std::string>(numEdgesSeen, std::string());

  abkwarn3(numEdgesSeen == numNets, "missing ", numNets - numEdgesSeen,
           " nets in .netD file");
  {
    int ii = -1;
    netD >> eathash(ii);
    abkwarn(!netD, "trailing lines in .netD file");
  }
}

void HGraphFixed::readAreM(const char* areMFileName) {
  ifstream areM(areMFileName);
  abkfatal2(areM, " Could not open ", areMFileName);

  string cName;      // cell name
  areM >> cName;     // skip cell name this time....
  areM >> cName;
  _numMultiWeights = 0;

  while (_numMultiWeights <= 100 && !cName.empty() && isdigit(cName[0])) {
    areM >> cName;
    _numMultiWeights++;
  }
  _numTotalWeights = _numMultiWeights;
  areM.close();
  areM.open(areMFileName);
  abkfatal(_numMultiWeights < 100, "error: there appear to be > 100 weights");
  abkfatal(_multiWeights.size() == 0, "multiWeights should be empty");

  _multiWeights = vector<HGWeight>(getNumNodes() * _numMultiWeights);

  int lineno = 0;

  while (areM >> cName) {
    if (cName.empty()) break;

    unsigned nodeIdx =
        parseNodeIndexToken(cName, _numTerminals, getNumNodes(), areMFileName);

    HGFNode& node = getNodeByIdx(nodeIdx);
    unsigned wtNum = 0;
    double nextWt;
    while (!areM.eof() && areM.peek() != '\n' && areM.peek() != '\r') {
      areM >> eatblank >> nextWt >> eatblank;
      setWeight(node.getIndex(), nextWt, wtNum++);
    }
    lineno++;
  }
  areM.close();
}

void HGraphFixed::readNodes(const char* nodesFileName) {
  ifstream nodes(nodesFileName);
  abkfatal2(nodes, " Could not open ", nodesFileName);
  cout << " Reading " << nodesFileName << " ... " << endl;

  int lineno = 1;
  nodes >> needcaseword("UCLA") >> needcaseword("nodes") >> needword("1.0") >>
      skiptoeol >> eathash(lineno);

  unsigned numNodes;
  nodes >> needcaseword("NumNodes") >> needword(":") >> my_isnumber(lineno) >>
      numNodes >> skiptoeol;
  lineno++;

  nodes >> needcaseword("NumTerminals") >> needword(":") >>
      my_isnumber(lineno) >> _numTerminals >> skiptoeol;
  lineno++;

  init(numNodes, 0);
  _nodeNamesMap.reserve(numNodes);

  string cName;
  vector<string> termNodeNames;
  vector<string> nonTermNodeNames;
  termNodeNames.reserve(_numTerminals);
  nonTermNodeNames.reserve(numNodes);

  nodes >> cName;
  for (unsigned nId = 0; nId < numNodes; nId++) {
    if (cName[0] == 0) {
      cerr << "Expected " << numNodes << " nodes, of them " << _numTerminals
           << " terminals " << endl;
      cerr << "Found only " << nId << " (through line " << lineno << ") "
           << endl;
      if (termNodeNames.size())
        cerr << "Last terminal name seen: " << termNodeNames.back() << ", ";
      if (nonTermNodeNames.size())
        cerr << "last non-terminal name seen: " << nonTermNodeNames.back();
      cerr << endl;
      abkfatal(0, "Not enough nodeNames in .nodes file");
    }

    string newNodeName(cName);

    bool isTerminal = false;

    nodes >> eatblank;

    while (!nodes.eof() && nodes.peek() != '\n' && nodes.peek() != '\r') {
      nodes >> cName;
      if (!strcasecmp(cName.c_str(), "terminal")) {
        isTerminal = true;
        nodes >> skiptoeol;
        break;
      } else
        nodes >> eatblank;
    }

    if (isTerminal)
      termNodeNames.push_back(newNodeName);
    else
      nonTermNodeNames.push_back(newNodeName);

    nodes >> skiptoeol;
    lineno++;
    nodes >> eathash(lineno);

    if (!nodes.eof())
      nodes >> cName;
    else
      cName[0] = 0;
  }

  if (_numTerminals != termNodeNames.size()) {
    cerr << "Expected " << _numTerminals << " terminals " << endl;
    cerr << "Found only " << termNodeNames.size() << " (through line " << lineno
         << ") " << endl;
    if (termNodeNames.size())
      cerr << "Last terminal name seen: " << termNodeNames.back() << ", ";
    if (nonTermNodeNames.size())
      cerr << "last non-terminal name seen: " << nonTermNodeNames.back();
    cerr << endl;
    abkfatal(0, "Did not find the expected number of terminals");
  }

  unsigned t;
  for (t = 0; t < termNodeNames.size(); t++) {
    string thisName = termNodeNames[t];
    _nodeNames.push_back(thisName);
    _nodeNamesMap[thisName] = t;
  }
  for (unsigned n = 0; n < nonTermNodeNames.size(); n++) {
    string thisName = nonTermNodeNames[n];
    _nodeNames.push_back(thisName);
    _nodeNamesMap[thisName] = t + n;
  }

  nodes.close();
}

void HGraphFixed::readNets(const char* netsFileName) {
  ifstream nets(netsFileName);
  abkfatal2(nets, " Could not open ", netsFileName);
  cout << " Reading " << netsFileName << " ... " << endl;

  int lineno = 1;
  nets >> needcaseword("UCLA") >> needcaseword("nets") >> needcaseword("1.0") >>
      skiptoeol;
  nets >> eathash(lineno);

  unsigned expectedNumPins = _numPins = 0;
  string cName;
  nets >> cName;
  if (cName == "NumNets") {
    nets >> skiptoeol;
    nets >> needcaseword("NumPins") >> needword(":") >> my_isnumber(lineno) >>
        expectedNumPins >> skiptoeol;
  } else if (cName == "NumPins")
    nets >> needword(":") >> my_isnumber(lineno) >> expectedNumPins >>
        skiptoeol;

#ifdef SIGNAL_DIRECTIONS
  _srcs.reserve(1 + (expectedNumPins / 3));
  _snks.reserve(1 + (expectedNumPins / 3));
  _srcSnks.reserve(1 + (expectedNumPins / 3));
#else
  _srcSnks.reserve(expectedNumPins);
#endif

  BitBoard nodesSeen(_nodes.size());
  string tmpNetName;

  unsigned netDegree;
  char dir;
  bool makingEdge;
  while (_numPins < expectedNumPins) {
    if (nets.eof()) {
      cerr << "Expected " << expectedNumPins << " got " << _numPins << endl;
      abkfatal(0, "mal-formed nets file (incorrect # pins)");
    }
    nets >> needcaseword("NetDegree") >> needword(":") >> my_isnumber(lineno) >>
        netDegree >> eatblank;
    char c = nets.peek();
    if (c != '\n' && c != '\r')
      nets >> tmpNetName >> skiptoeol;
    else
      tmpNetName.clear();
    HGFEdge* edge;
    if (netDegree > 1) {
      edge = addEdge();
      if (!tmpNetName.empty()) {
        _netNames.push_back(tmpNetName);
        _netNamesMap[tmpNetName] = edge->getIndex();
      } else
        _netNames.push_back(string());
      abkassert(edge->getWeight() > 0, "constructed a 0-weight edge");
      makingEdge = true;
    } else {
      cerr << " Net of degree 1 found ";
      if (!tmpNetName.empty()) cerr << " named " << tmpNetName;
      abkfatal(0, "This parser does not allow nets of degree < 2");
      makingEdge = false;
    }

    nodesSeen.clear();

    for (unsigned n = 0; n < netDegree; n++) {
      if (!makingEdge) {
        nets >> skiptoeol;
        _numPins++;
        continue;
      }

      nets >> cName;
      abkfatal(cName != "NetDegree",
               "did not find expected number of node names");
      abkfatal2(haveSuchNode(cName.c_str()), "no such node: ", cName.c_str());
      HGFNode& node = getNodeByName(cName.c_str());
      unsigned nId = node.getIndex();

      if (nodesSeen.isBitSet(nId)) {
        _numPins++;
        nets >> skiptoeol;
        nets >> eathash(lineno);
        lineno++;
        continue;
      }

      nets >> dir >> skiptoeol;
      dir = toupper(dir);
      switch (dir) {
        case ('I'):
#ifdef SIGNAL_DIRECTIONS
          addSnk(node, *edge);
          break;
#endif
        case ('O'):
#ifdef SIGNAL_DIRECTIONS
          addSrc(node, *edge);
          break;
#endif
        case ('B'):
          addSrcSnk(node, *edge);
          break;
        default:  // error
          abkfatal3(0, dir, ": unexpected symbol in line ", lineno);
      };
      _numPins++;
      lineno++;
      nets >> eathash(lineno);
      nodesSeen.setBit(nId);

      if (_param.verb.getForActions() > 10)
        cout << "Parsed net " << edge->getIndex() << " degree " << netDegree
             << endl;
    }
    nets >> eathash(lineno);
  }

  nets.close();
}

void HGraphFixed::readWts(const char* wtsFileName) {
  ifstream wts(wtsFileName);
  abkfatal2(wts, " Could not open ", wtsFileName);
  cout << " Reading " << wtsFileName << " ... " << endl;

  int lineno = 1;

  wts >> needcaseword("UCLA") >> needcaseword("wts") >> needword("1.0") >>
      skiptoeol;
  wts >> eathash(lineno);

  _numMultiWeights = 0;
  string cName;      // cell name

  if (!wts.eof()) {
    wts >> cName;  // skip cell name this time....
    if (!wts.eof()) {
      wts >> cName;
      while (!wts.eof() && _numMultiWeights <= 100 && !cName.empty() &&
             isdigit(cName[0])) {
        wts >> cName;
        _numMultiWeights++;
      }
    }
  }

  _numTotalWeights = _numMultiWeights;
  wts.close();
  abkfatal(_numMultiWeights < 100,
           "error: there appear to be > 100 weights OR cell names start with "
           "numbers");
  abkfatal(_multiWeights.size() == 0, "multiWeights should be empty");

  if (_numTotalWeights == 0) {  // set unit weights
    loadUnitWeights();
    return;
  }

  _multiWeights = vector<HGWeight>(getNumNodes() * _numMultiWeights);

  lineno = 0;

  ifstream wts2(wtsFileName);
  abkfatal2(wts2, " Could not open ", wtsFileName);

  wts2 >> needcaseword("UCLA") >> needcaseword("wts") >> needword("1.0") >>
      skiptoeol;
  wts2 >> eathash(lineno);

  int NodeNetNameWarningCount = 0;

  while (!wts2.eof() && wts2 >> eathash(lineno) >> cName) {
    if (!wts2.eof()) {
      if (cName.empty()) break;
      if (haveSuchNode(cName.c_str()) && haveSuchNet(cName.c_str()) &&
          NodeNetNameWarningCount < 10) {
        NodeNetNameWarningCount++;
        string msg = "Node and Net with the same name \"" + cName +
                     "\",\n\t .wts preference given to Net.\n";
        if (NodeNetNameWarningCount == 10)
          msg +=
              "This warning has occurred 10 times, and will not be repeated.\n";
        abkwarn(0, msg.c_str());
      }

      if (haveSuchNet(cName.c_str())) {
        HGFEdge& net = getNetByName(cName.c_str());
        double nextWt;
        wts2 >> eatblank >> nextWt >> eatblank;
        net.setWeight(nextWt);
      } else if (haveSuchNode(cName.c_str())) {
        HGFNode& node = getNodeByName(cName.c_str());
        unsigned wtNum = 0;
        double nextWt;
        while (!wts2.eof() && wts2.peek() != '\n' && wts2.peek() != '\r') {
          wts2 >> eatblank >> nextWt >> eatblank;
          setWeight(node.getIndex(), nextWt, wtNum++);
        }
      }
      lineno++;
      cName.clear();
    }
  }

  wts2.close();
}
