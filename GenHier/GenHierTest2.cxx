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

// Interactive hierarchy construction driver.
// It reads a delimiter set and a list of names from stdin, builds the
// hierarchy, and writes the result to test2.hcl.

#include <iostream>
#include <string>
#include <vector>

#include "genHier.h"

using std::cin;
using std::cout;
using std::endl;
using std::string;
using std::vector;

int main() {
  cout << "Enter a string of one-character delimiters (e.g., '/|_'):  ";
  string delim;
  cin >> delim;
  cout << "The delimiter set is " << delim << endl;
  cout << "Now enter a list of words containing those delimiters "
       << "\n (first empty line ends input) " << endl;
  vector<string> nameStorage;
  string word;
  while (cin >> word) {
    if (word.empty()) break;
    nameStorage.push_back(word);
  }
  vector<const char*> names;
  names.reserve(nameStorage.size());
  for (unsigned i = 0; i < nameStorage.size(); ++i) {
    names.push_back(nameStorage[i].c_str());
  }
  if (names.size() > 1)
    cout << "Read " << names.size() << " words, will now ctruct a hierarhy...";
  else {
    cout << "Need two or more words" << endl;
    return 1;
  }
  GenericHierarchy gHie(names, delim.c_str());
  cout << " done \n saving test2.hcl ...";
  gHie.saveHCL("test2.hcl");
  cout << " done " << endl;
}
