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

#ifdef _MSC_VER
#pragma warning(disable : 4786)
#endif

#include <ABKCommon/abkio.h>
#include <ABKCommon/infolines.h>

#include <Placement/placeWOri.h>

#include <errno.h>
#include <limits.h>
#include <stdlib.h>
#include <string>

using std::endl;
using std::ifstream;
using std::ofstream;
using std::ostream;
using std::string;
using std::vector;
using std::setw;

namespace {

static void initOrients(Orientation* orients, unsigned nPts,
                        const Orientation& ori) {
  if (!orients) return;
  for (unsigned i = 0; i < nPts; ++i) orients[i] = ori;
}

static void copyOrients(Orientation* dest, const Orientation* src,
                        unsigned nPts) {
  if (!dest || !src) return;
  for (unsigned i = 0; i < nPts; ++i) dest[i] = src[i];
}

static unsigned parseNodeNameSuffix(const char* str) {
  abkfatal(str != NULL && str[0] != '\0' && str[1] != '\0',
           "Malformed node name");
  errno = 0;
  char* end = NULL;
  const unsigned long parsed = strtoul(str + 1, &end, 10);
  abkfatal(end != str + 1 && *end == '\0' && errno != ERANGE &&
               parsed <= UINT_MAX,
           "Malformed node name index");
  return static_cast<unsigned>(parsed);
}

}  // namespace

static unsigned getNodeIndexFromName(const char* str, unsigned numTerms) {
  const unsigned indexInName = parseNodeNameSuffix(str);
  const char prefix = toupper(static_cast<unsigned char>(str[0]));
  if (prefix == 'P') {
    abkfatal(indexInName > 0, "Terminal node names are one-based");
    const unsigned index = indexInName - 1;
    abkfatal2(index < numTerms, index + 1, ": terminal index too big");
    return index;
  } else if (prefix == 'A') {
    abkfatal(indexInName <= UINT_MAX - numTerms, "Node index is too large");
    return indexInName + numTerms;
  } else
    abkfatal2(0, str[0], " encountered instead of 'a' or 'p'");
  return 0;
}

static void getNameFromNodeIndex(unsigned index, unsigned numTerms, char* str) {
  if (index < numTerms) {
    str[0] = 'p';
    ++index;
  } else {
    str[0] = 'a';
    index -= numTerms;
  }
  snprintf(str + 1, 29, "%u", index);
}

PlacementWOrient::PlacementWOrient(unsigned n) : Placement(n), orients(NULL) {
  if (nPts) orients = new Orientation[n];
}

PlacementWOrient::PlacementWOrient(unsigned n, Point pt, Orientation ori)
    : Placement(n, pt) {
  if (nPts) {
    orients = new Orientation[nPts];
    initOrients(orients, nPts, ori);
  } else
    orients = NULL;
}

PlacementWOrient::PlacementWOrient(const PlacementWOrient& pl) : Placement(pl) {
  if (nPts) {
    orients = new Orientation[nPts];
    copyOrients(orients, pl.orients, nPts);
  } else
    orients = NULL;
}

PlacementWOrient::PlacementWOrient(const vector<Point>& pts)
    : Placement(pts), orients(NULL) {
  if (nPts) {
    orients = new Orientation[nPts];
    initOrients(orients, nPts, Orientation());
  }
}

PlacementWOrient::PlacementWOrient(const vector<Point>& pts,
                                   const vector<Orientation>& ori)
    : Placement(pts), orients(NULL) {
  abkfatal(pts.size() == ori.size(), "points and orients must have same size");
  if (nPts) orients = new Orientation[nPts];
  copyOrients(orients, ori.size() ? &ori[0] : 0, nPts);
}

PlacementWOrient::PlacementWOrient(const Placement& pl,
                                   const vector<Orientation>& ori)
    : Placement(pl), orients(NULL) {
  abkfatal(ori.size() == nPts, "points and orients must have same size");
  if (nPts) orients = new Orientation[nPts];
  copyOrients(orients, ori.size() ? &ori[0] : 0, nPts);
}

PlacementWOrient::PlacementWOrient(const Mapping& pullBackMap,
                                   const PlacementWOrient& from)
    : Placement(pullBackMap.getSourceSize()), orients(NULL) {
  if (nPts == 0) {
    points = NULL;
    return;
  }
  points = new Point[nPts];
  orients = new Orientation[nPts];
  for (unsigned k = 0; k < nPts; k++) {
    points[k] = from[pullBackMap[k]];
    orients[k] = from.getOrient(pullBackMap[k]);
  }
}

PlacementWOrient::PlacementWOrient(const char* plFileName, unsigned nTerms,
                                   unsigned nCore)
    : Placement(nTerms + nCore) {
  abkfatal(nPts > 0, "can't construct an empty PlaceWOrient from a file");
  orients = new Orientation[nPts];

  ifstream in(plFileName);
  int lineno = 1;

  in >> needcaseword("UCLA") >> needcaseword("pl") >> needword("1.0") >>
      skiptoeol;
  in >> eathash(lineno);

  string cName;
  Point pt;
  string oriStr;
  while (!in.eof()) {
    in >> cName >> pt.x >> pt.y >> eatblank;
    unsigned nodeId = getNodeIndexFromName(cName.c_str(), nTerms);
    abkfatal(nodeId < nPts, "Node exceeds nPts");
    points[nodeId] = pt;

    if (in.peek() == ':') {
      in >> needword(":") >> oriStr;
      orients[nodeId] = Orientation(oriStr.c_str());
    }

    in >> skiptoeol;
  }
  in.close();
}

void PlacementWOrient::reorder(const Permutation& pm) {
  abkassert(pm.getSize() == getSize(),
            "Can't reorder placement: wrong size of permutation\n");
  Point* newPoints = new Point[getSize()];
  Orientation* newOris = new Orientation[getSize()];
  unsigned k;
  for (k = 0; k < getSize(); k++) {
    newPoints[pm[k]] = points[k];
    newOris[pm[k]] = orients[k];
  }
  copyOrients(orients, newOris, getSize());
  for (k = 0; k < getSize(); k++) points[k] = newPoints[k];
  delete[] newPoints;
  delete[] newOris;
}
void PlacementWOrient::reorderBack(const Permutation& pm) {
  abkassert(pm.getSize() == getSize(),
            "Can't reorder placement: wrong size of permutation\n");
  Point* newPoints = new Point[getSize()];
  Orientation* newOris = new Orientation[getSize()];

  unsigned k;
  for (k = 0; k < getSize(); k++) {
    newPoints[k] = points[pm[k]];
    newOris[k] = orients[pm[k]];
  }
  for (k = 0; k < getSize(); k++) points[k] = newPoints[k];
  copyOrients(orients, newOris, getSize());
  delete[] newPoints;
  delete[] newOris;
}

PlacementWOrient& PlacementWOrient::operator=(const PlacementWOrient& from) {
  if (from.points == points && from.orients == orients) return *this;

  delete[] points;
  delete[] orients;
  nPts = from.size();
  if (nPts) {
    points = new Point[nPts];
    orients = new Orientation[nPts];
    for (unsigned i = 0; i < getSize(); i++) {
      points[i] = from[i];
      orients[i] = from.getOrient(i);
    }
  } else {
    points = NULL;
    orients = NULL;
  }
  return *this;
}

void PlacementWOrient::save(const char* plFileName, unsigned numTerms) const {
  ofstream out(plFileName);

  out << "UCLA pl 1.0" << endl;
  out << TimeStamp() << User() << Platform() << endl << endl;

  char cName[30];
  for (unsigned i = 0; i < nPts; i++) {
    getNameFromNodeIndex(i, numTerms, cName);
    out << setw(8) << cName << "  " << setw(10) << points[i].x << "  "
        << setw(10) << points[i].y << " : " << orients[i] << endl;
  }

  out.close();
}

ostream& operator<<(ostream& out, const PlacementWOrient& arg) {
  TimeStamp tm;
  out << tm;
  out << "Total points : " << arg.nPts << endl;
  for (unsigned i = 0; i < arg.nPts; i++) out << arg.points[i] << " " << "\n";
  return out;
}
