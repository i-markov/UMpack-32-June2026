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

#include <algorithm>
#include <string>
#include <unordered_map>
#include <unordered_set>

#include "ABKCommon/abk_hash_common.h"
#include "genHier.h"

using std::min;
using std::sort;
using std::string;
using std::vector;

struct GenHier_StringHash {
  size_t operator()(const string& str) const {
    return ABKCStringHash()(str.c_str());
  }
};

typedef std::unordered_set<string, GenHier_StringHash> GenHier_StringSet;

void gen_hier_common_prefix(const char* a, const char* b, const char* delim_set,
                            char* commonPrefix) {
  abkfatal(a != NULL && b != NULL,
           "gen_hier_common_prefix requires non-null input strings");
  abkfatal(commonPrefix != NULL,
           "gen_hier_common_prefix requires caller-owned output storage");

  if (a[0] == '\0' || b[0] == '\0') {
    strcpy(commonPrefix, "");
    return;
  }
  unsigned lena = strlen(a), lenb = strlen(b), lenx = 0;
  if (delim_set == NULL || delim_set[0] == '\0') {
    for (; lenx != min(lena, lenb); lenx++) {
      if (a[lenx] != b[lenx]) break;
    }

    strncpy(commonPrefix, a, lenx);
    commonPrefix[lenx] = '\0';
    return;
  }
  const char *oa = strpbrk(a, delim_set), *ob = strpbrk(b, delim_set);
  lena = oa ? oa - a : strlen(a);
  lenb = ob ? ob - b : strlen(b);
  if (lena != lenb || strncmp(a, b, lena)) {
    strcpy(commonPrefix, "");
    return;
  } else if (!oa || !ob) {
    strncpy(commonPrefix, a, lena);
    commonPrefix[lena] = '\0';
    return;
  }

  do {
    lenx = lena;
    oa = strpbrk(oa + 1, delim_set);
    ob = strpbrk(ob + 1, delim_set);
    lena = oa ? oa - a : strlen(a);
    lenb = ob ? ob - b : strlen(b);
    if (!oa || !ob) break;
  } while (lena == lenb && strncmp(a, b, lena) == 0);

  if (oa == NULL) lena = strlen(a);
  if (ob == NULL) lenb = strlen(b);
  if (lena == lenb && strncmp(a, b, lena) == 0) lenx = lena;

  strncpy(commonPrefix, a, lenx);
  commonPrefix[lenx] = '\0';

  return;
}

static string genHierCommonPrefixString(const string& a, const string& b,
                                        const char* delimiter_set) {
  const unsigned maxLen = a.size() > b.size() ? a.size() : b.size();
  vector<char> commonPrefix(maxLen + 1, '\0');
  gen_hier_common_prefix(a.c_str(), b.c_str(), delimiter_set,
                         &commonPrefix[0]);
  return string(&commonPrefix[0]);
}

GenericHierarchy::GenericHierarchy(const vector<const char*>& nodeNames,
                                   const char* delimiter_set)
    : _names(nodeNames.size()), _root(UINT_MAX) {
  vector<string> tmpNames;
  tmpNames.reserve(nodeNames.size());
  unsigned numLeaves = nodeNames.size();
  GenHier_StringSet leaves;
  unsigned i;
  for (i = 0; i < numLeaves; i++) {
    tmpNames.push_back(nodeNames[i]);
    leaves.insert(tmpNames.back());
  }
  sort(tmpNames.begin(), tmpNames.end());

  GenHier_StringSet addedClusters;
  for (i = 1; i != numLeaves; i++) {
    const string prefix =
        genHierCommonPrefixString(tmpNames[i - 1], tmpNames[i], delimiter_set);
    if (leaves.find(prefix) != leaves.end()) continue;
    const string clusterName = prefix.empty() ? "ROOT" : prefix;
    if (addedClusters.find(clusterName) == addedClusters.end()) {
      addedClusters.insert(clusterName);
      tmpNames.push_back(prefix);
    }
  }
  sort(tmpNames.begin(), tmpNames.end());

  std::unordered_map<string, string, GenHier_StringHash> parentNames;
  GenHier_StringSet missingChildren;

  if (tmpNames[0].empty()) {
    parentNames["ROOT"] = "ROOT";
  } else
    parentNames[tmpNames[0]] = tmpNames[0];

  unsigned nameSize = tmpNames.size();
  for (i = 1; i != nameSize; i++) {
    const string prefix =
        genHierCommonPrefixString(tmpNames[i - 1], tmpNames[i], delimiter_set);
    if (leaves.find(prefix) == leaves.end()) {
      parentNames[tmpNames[i]] = prefix.empty() ? "ROOT" : prefix;
    } else {
      string clusterName(prefix);
      clusterName += "++cluster";
      parentNames[tmpNames[i]] = clusterName;
      missingChildren.insert(prefix);
    }
  }

  for (i = 0; i < numLeaves; i++) {
    _names[i] = nodeNames[i];
    _namesToId[_names[i]] = i;
  }

  GenHier_StringSet::iterator it;

  for (it = addedClusters.begin(); it != addedClusters.end(); it++) {
    _names.push_back(*it);
    _namesToId[_names.back()] = i++;
  }

  for (it = missingChildren.begin(); it != missingChildren.end(); it++) {
    string clusterName(*it);
    clusterName += "++cluster";
    parentNames[clusterName] = parentNames[*it];
    parentNames[*it] = clusterName;
    _names.push_back(clusterName);
    _namesToId[_names.back()] = i++;
  }

  // from the parentNames, generate the parent and then children Id's
  _parents = vector<unsigned>(_names.size(), GENH_DELETED_NODE);
  _children = vector<vector<unsigned> >(_names.size());

  for (i = 0; i < _names.size(); i++) {
    const string parentName = parentNames[_names[i]];
    if (_namesToId.count(parentName) == 0) {
      if (_names[i] != parentName)
        abkfatal3(false, _names[i].c_str(),
                  "'s parent was not defined: ", parentName.c_str());
      _parents[i] = i;
      _root = i;
      continue;
    }

    _parents[i] = _namesToId[parentName];

    if (_parents[i] == i)
      _root = i;
    else
      _children[_parents[i]].push_back(i);
  }

  /*  for(i=0; i!=_names.size(); i++)
      {
         cout <<"\"" << _names[i] << "\"   \""
                     << parentNames[_names[i]]<< "\"" << endl;
      }
  */
}
