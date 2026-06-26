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

// 980217  dv added orientation-SET for cells.
#ifdef _MSC_VER
#pragma warning(disable : 4786)
#endif

#include <Placement/symOri.h>

#include <cstring>
#include <iomanip>
#include <iostream>

using std::cout;
using std::endl;
using std::ostream;

namespace {

const char* kSymmetryText[8] = {" ",     " X ",     " Y ",     " Y X ",
                                " R90 ", " R90 X ", " R90 Y ", " R90 Y X "};

const char* kOrientText[8] = {"FN ", " N ", "FE ", " E ",
                              "FS ", " S ", "FW ", " W "};

const char* kOrientNames[8] = {"N", "FN", "E", "FE", "S", "FS", "W", "FW"};

}  // namespace

Symmetry::Symmetry(const char* txt) {
  if (strstr(txt, "R90") || strstr(txt, "r90"))
    rot90 = true;
  else
    rot90 = false;
  if (strchr(txt, 'Y') || strchr(txt, 'y'))
    y = true;
  else
    y = false;
  if (strchr(txt, 'X') || strchr(txt, 'x'))
    x = true;
  else
    x = false;
}

Symmetry::operator const char*() const {
  return kSymmetryText[unsigned(*this)];
}

Orient::Orient(const char* txt) {
  _angle = 0;  // initialize to N orientation
  if (strchr(txt, 'F') || strchr(txt, 'f'))
    _faceUp = false;
  else
    _faceUp = true;
  if (strchr(txt, 'N') || strchr(txt, 'n')) {
    _angle = 0;
    return;
  }
  if (strchr(txt, 'E') || strchr(txt, 'e')) {
    _angle = 1;
    return;
  }
  if (strchr(txt, 'S') || strchr(txt, 's')) {
    _angle = 2;
    return;
  }
  if (strchr(txt, 'W') || strchr(txt, 'w')) {
    _angle = 3;
    return;
  }
}

Orient::operator const char*() const {
  return kOrientText[_angle * 2 + (_faceUp ? 1 : 0)];
}

void OrientationSet::setBits(unsigned x) {
  N = x % 2;
  x /= 2;
  FN = x % 2;
  x /= 2;
  E = x % 2;
  x /= 2;
  FE = x % 2;
  x /= 2;
  S = x % 2;
  x /= 2;
  FS = x % 2;
  x /= 2;
  W = x % 2;
  x /= 2;
  FW = x % 2;
}

unsigned OrientationSet::bits() const {
  return N + 2 * FN + 4 * E + 8 * FE + 16 * S + 32 * FS + 64 * W + 128 * FW;
}

Orient OrientationSet::orientFromIndex(unsigned idx) {
  abkassert(idx < 8, "Invalid orientation index");
  return Orient(kOrientNames[idx]);
}

OrientationSet::OrientationSet(const Symmetry& sym) {
  *this = OrientationSet(sym, Orient("N"));
}

OrientationSet::OrientationSet(const Symmetry& sym, const Orient& orient)
/*
   Convert the symmetry, orientation pair into an orientation set.

   It would be nice to do something slick here like:
       Start with the single orientation, then generate the
       whole set by applying symmetry generation operations...

   But I don't think there's an easy way to reverse the bits
   in an unsigned, and Symmetry already exposes its
   internal representation, so we'll just look at the bits
   and do a big switch.  Sigh.
*/
{
  unsigned r;
  unsigned ori = orient.getAngle() + (orient.isFaceUp() ? 1 : 0);

  //  cout << "sym: " << sym << " orient " << orient << endl;
  //  cout << "sym, or: " << unsigned(sym) << " " << or << endl;

  switch (unsigned(sym)) {
    case 0:  // No symmetry
      switch (ori) {
        case 1:
          r = 1;  // N
          break;
        case 0:
          r = 2;  // FN
          break;
        case 91:
          r = 4;  // E
          break;
        case 90:
          r = 8;  // FE
          break;
        case 181:
          r = 16;  // S
          break;
        case 180:
          r = 32;  // FS
          break;
        case 271:
          r = 64;  // W
          break;
        case 270:
          r = 128;  // FW
          break;
        default:
          cout << "bad value case 0: " << ori << endl;
      }
    case 1:  // X only
      switch (ori) {
        case 1:
        case 180:
          r = 33;  // N+FS
          break;
        case 0:
        case 181:
          r = 18;  // FN+S
          break;
        case 90:
        case 91:
          r = 12;  // E+FE
          break;
        case 270:
        case 271:
          r = 192;  // W+FW
          break;
        default:
          cout << "bad value case 1: " << ori << endl;
      }
      break;
    case 2:  // Y only
      switch (ori) {
        case 0:
        case 1:
          r = 3;  // N+FN
          break;
        case 91:
        case 270:
          r = 132;  // E+FW
          break;
        case 90:
        case 271:
          r = 72;  // FE+W
          break;
        case 180:
        case 181:
          r = 48;  // S+FS
          break;
        default:
          cout << "bad value case 2: " << ori << endl;
      }
      break;
    case 3:  // X and Y
      switch (ori) {
        case 0:
        case 1:
        case 180:
        case 181:
          r = 51;  // N+FN+S+FS
          break;
        default:
          r = 214;  // E+FE+W+FW
      }
      break;
    case 4:  // Rot90 only
      switch (ori) {
        case 1:
        case 91:
        case 181:
        case 271:
          r = 85;  // N+E+S+W
          break;
        default:
          r = 170;  // FN+FE+FS+FW
      }
      break;

    default:  // all other cases are the same
      r = 255;
  }

  // Now just set all the bits

  setBits(r);
}

OrientationSet::OrientationSet(unsigned x) { setBits(x); }

Orient OrientationSet::getNth(unsigned n) {
  const unsigned mask = bits();
  for (unsigned idx = 0; idx < 8; ++idx)
    if (mask & (1u << idx)) {
      if (n == 0) return orientFromIndex(idx);
      --n;
    }
  abkfatal(false, "Requested non-existent member of orientation set.");

  return Orient(NULL);  // satisfy compiler
}

Orient OrientationSet::getOrient() {
  abkassert(isNotEmpty(), "Tried to extract orientation from empty set");
  return getNth(0);

  return Orient(NULL);  // satisfy compiler
}

unsigned OrientationSet::size() { return N + FN + E + FE + S + FS + W + FW; }

ostream& operator<<(ostream& out, const OrientationSet& os) {
  out << "OrientationSet{";
  if (os.N) out << "N ";
  if (os.FN) out << "FN ";
  if (os.E) out << "E ";
  if (os.FE) out << "FE ";
  if (os.S) out << "S ";
  if (os.FS) out << "FS ";
  if (os.W) out << "W ";
  if (os.FW) out << "FW ";
  out << "}";
  return out;
}

ostream& operator<<(ostream& out, const Symmetry& sym) {
  out << (const char*)sym;
  return out;
}

ostream& operator<<(ostream& out, const Orient& ori) {
  out << (const char*)ori;
  return out;
}
