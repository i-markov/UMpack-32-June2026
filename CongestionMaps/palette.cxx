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
#include "palette.h"

#include <iomanip>
#include <iostream>
#include <sstream>
#include <utility>

#include "ABKCommon/abkassert.h"

using std::cout;
using std::endl;
using std::hex;
using std::make_pair;
using std::string;
using std::vector;
using std::stringstream;

PaletteBox* PaletteBox::_impl = 0;

Palette PaletteBox::Hot(void) {
  return make_pair(getImpl()->_hotMap, getImpl()->_hotColors);
}

Palette PaletteBox::GoodAndBad(void) {
  return make_pair(getImpl()->_goodAndBadMap, getImpl()->_goodAndBadColors);
}

Palette PaletteBox::Greyscale(void) {
  return make_pair(getImpl()->_greyscaleMap, getImpl()->_greyscaleColors);
}

PaletteBox* PaletteBox::getImpl(void) {
  if (!_impl) _impl = new PaletteBox();
  return _impl;
}

PaletteBox::PaletteBox(void) {
  // some local variables with hard coded initializations
  // to initialize vectors and strings with hardcoded values

  const char* Colors[] = {"a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k",
                          "l", "m", "n", "o", "p", "q", "r", "s", "t", "u", "v",
                          "w", "x", "y", "z", "A", "B", "C", "D", "E", "F", "G",
                          "H", "I", "J", "K", "L", "M", "N", "O", "P", "Q", "R",
                          "S", "T", "U", "V", "W", "X", "Y", "Z", "0", "1", "2",
                          "3", "4", "5", "6", "7", "8", "9", ",", "."};

  /* matlab code to generate a proper palette
    newpalette = round(sqrt(hot(64))*255);
    i = 1;
    colorlist = '';
    while(i <= 64)
      j = 1;
      color = '';
      while(j <= 3)
        color = strcat(color,sprintf('%02x',newpalette(i,j)));
        j = j + 1;
      end
      colorlist = strvcat(colorlist,color);
      i = i + 1;
    end
  */

  const char* hotMap[] = {
      /* columns rows colors chars-per-pixel */
      "a c #340000", "b c #4a0000", "c c #5a0000", "d c #680000", "e c #740000",
      "f c #800000", "g c #8a0000", "h c #930000", "i c #9c0000", "j c #a50000",
      "k c #ad0000", "l c #b40000", "m c #bc0000", "n c #c30000", "o c #ca0000",
      "p c #d00000", "q c #d70000", "r c #dd0000", "s c #e30000", "t c #e90000",
      "u c #ef0000", "v c #f40000", "w c #fa0000", "x c #ff0000", "y c #ff3400",
      "z c #ff4a00", "A c #ff5a00", "B c #ff6800", "C c #ff7400", "D c #ff8000",
      "E c #ff8a00", "F c #ff9300", "G c #ff9c00", "H c #ffa500", "I c #ffad00",
      "J c #ffb400", "K c #ffbc00", "L c #ffc300", "M c #ffca00", "N c #ffd000",
      "O c #ffd700", "P c #ffdd00", "Q c #ffe300", "R c #ffe900", "S c #ffef00",
      "T c #fff400", "U c #fffa00", "V c #ffff00", "W c #ffff40", "X c #ffff5a",
      "Y c #ffff6e", "Z c #ffff80", "0 c #ffff8f", "1 c #ffff9c", "2 c #ffffa9",
      "3 c #ffffb4", "4 c #ffffbf", "5 c #ffffca", "6 c #ffffd3", "7 c #ffffdd",
      "8 c #ffffe6", "9 c #ffffef", ", c #fffff7", ". c #ffffff"};

  // compute the goodAndBad
  vector<string> goodAndBadMap;
  for (unsigned i = 0; i < NUM_COLORS / 2; ++i) {
    unsigned red = 0, green = 0, blue = 0;
    green = 255 * (NUM_COLORS / 2 - i) / (NUM_COLORS / 2);
    stringstream ss;
    ss << Colors[i] << " c #" << get2CharHexNum(red) << get2CharHexNum(green)
       << get2CharHexNum(blue);
    goodAndBadMap.push_back(ss.str().c_str());
  }
  for (unsigned i = 0; i < NUM_COLORS / 2; ++i) {
    unsigned red = 0, green = 0, blue = 0;
    red = 255 * (i) / (NUM_COLORS / 2);
    stringstream ss;
    ss << Colors[i + NUM_COLORS / 2] << " c #" << get2CharHexNum(red)
       << get2CharHexNum(green) << get2CharHexNum(blue);
    goodAndBadMap.push_back(ss.str().c_str());
  }

  // compute the greyscale
  vector<string> greyscaleMap;
  for (unsigned i = 0; i < NUM_COLORS; ++i) {
    unsigned red = 0, green = 0, blue = 0;
    red = green = blue = 255 * (i) / (NUM_COLORS);
    stringstream ss;
    ss << Colors[i] << " c #" << get2CharHexNum(red) << get2CharHexNum(green)
       << get2CharHexNum(blue);
    greyscaleMap.push_back(ss.str().c_str());
  }

  for (unsigned i = 0; i < NUM_COLORS; ++i) _hotMap.push_back(hotMap[i]);
  for (unsigned i = 0; i < NUM_COLORS; ++i) _hotColors.push_back(Colors[i]);

  for (unsigned i = 0; i < NUM_COLORS; ++i)
    _goodAndBadMap.push_back(goodAndBadMap[i]);
  for (unsigned i = 0; i < NUM_COLORS; ++i)
    _goodAndBadColors.push_back(Colors[i]);

  for (unsigned i = 0; i < NUM_COLORS; ++i)
    _greyscaleMap.push_back(greyscaleMap[i]);
  for (unsigned i = 0; i < NUM_COLORS; ++i)
    _greyscaleColors.push_back(Colors[i]);
}

string PaletteBox::get2CharHexNum(unsigned num) {
  stringstream ss;
  ss << hex << num;
  string rval = ss.str().c_str();
  if (rval.size() == 1) rval = "0" + rval;
  abkassert(rval.size() == 2, "Must return 2 chars");
  return rval;
}
