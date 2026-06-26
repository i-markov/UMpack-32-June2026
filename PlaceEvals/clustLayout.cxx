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

#include "ABKCommon/abkassert.h"
#include "ABKCommon/abkio.h"
#include "clustLayout.h"

using std::endl;
using std::ifstream;
using std::ostream;
using std::vector;

void ClusteredLayout::print(ostream &o) {
  for (unsigned i = 0; i < _numLevels; i++) {
    o << "Level " << i << endl;
    for (unsigned j = 0; j < _layoutTree[i].size(); j++) {
      o << "  Rect " << j << " : " << _layoutTree[i][j];
      o << "  Parent : " << _parentMap[i][j] << endl;
      o << "  Children : " << _childMap[i][j].size() << endl;
      for (unsigned k = 0; k < _childMap[i][j].size(); k++)
        o << "        " << _childMap[i][j][k] << endl;
    }
  }
}

ClusteredLayout::ClusteredLayout(const char *fileList) {
  ifstream inFile(fileList);
  abkfatal(inFile, "cannot open file for input");
  int lineNum = 1;
  unsigned i;
  inFile >> eathash(lineNum) >> my_isnumber(lineNum) >> _numLevels >> skiptoeol;
  for (i = 0; i < _numLevels; i++) {
    char fileName[200];
    inFile >> eathash(lineNum) >> fileName >> skiptoeol;
    LayoutBBoxes bboxList(fileName);
    _layoutTree.push_back(bboxList);

    vector<vector<unsigned> > children(bboxList.size());
    _childMap.push_back(children);

    vector<unsigned> parent(bboxList.size(), bboxList.size());
    _parentMap.push_back(parent);
  }
  inFile.close();

  for (unsigned level = 1; level < _numLevels; level++) {
    vector<vector<unsigned> > &children = _childMap[level - 1];
    vector<unsigned> &parent = _parentMap[level];
    LayoutBBoxes &currentLevel = _layoutTree[level],
                 &previousLevel = _layoutTree[level - 1];

    for (i = 0; i < currentLevel.size(); i++)
      for (unsigned j = 0; j < previousLevel.size(); j++) {
        if (previousLevel[j].contains(currentLevel[i])) {
          children[j].push_back(i);
          parent[i] = j;
        }
      }
  }
}
