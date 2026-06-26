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

// created by Stefanus Mantik on 06/2/98
#ifdef _MSC_VER
#pragma warning(disable : 4786)
#endif

#include "expectedBBox.h"

#include "Geoms/bbox.h"

using std::swap;

static void HBBoxShortToLeft(double& cLow, double& cHigh, double nLow,
                             double nHigh) {
  if (nLow < cLow) {
    swap(cLow, nLow);
    swap(cHigh, nHigh);
  }
  double cLen = cHigh - cLow, nLen = nHigh - nLow;
  if (nLow >= cHigh) return;
  double distance = cHigh - nLow;
  if (nHigh > cHigh)
    cHigh -= distance * distance * distance / (3 * cLen * nLen);
  else
    cHigh -=
        ((nLen * nLen / 3) - (nLen * distance) + (distance * distance)) / cLen;
  if (cHigh < cLow) cHigh = cLow;
}

static void HBBoxShortToRight(double& cLow, double& cHigh, double nLow,
                              double nHigh) {
  if (nHigh > cHigh) {
    swap(cLow, nLow);
    swap(cHigh, nHigh);
  }
  double cLen = cHigh - cLow, nLen = nHigh - nLow;
  if (nHigh <= cLow) return;
  double distance = nHigh - cLow;
  if (nLow < cLow)
    cLow += distance * distance * distance / (3 * cLen * nLen);
  else
    cLow +=
        ((nLen * nLen / 3) - (nLen * distance) + (distance * distance)) / cLen;
  if (cHigh < cLow) cLow = cHigh;
}

BBox HBBoxUpdateTopLeft(const BBox& current, const BBox& next) {
  BBox newBox;
  // set the left coordinate
  double cLow = current.xMin, cHigh = current.xMax;
  double nLow = next.xMin, nHigh = next.xMax;
  HBBoxShortToLeft(cLow, cHigh, nLow, nHigh);
  newBox.xMin = cLow;
  newBox.xMax = cHigh;

  // set the top coordinate
  cLow = current.yMin;
  cHigh = current.yMax;
  nLow = next.yMin;
  nHigh = next.yMax;
  HBBoxShortToRight(cLow, cHigh, nLow, nHigh);
  newBox.yMin = cLow;
  newBox.yMax = cHigh;
  return newBox;
}

BBox HBBoxUpdateBottomRight(const BBox& current, const BBox& next) {
  BBox newBox;
  // set the right coordinate
  double cLow = current.xMin, cHigh = current.xMax;
  double nLow = next.xMin, nHigh = next.xMax;
  HBBoxShortToRight(cLow, cHigh, nLow, nHigh);
  newBox.xMin = cLow;
  newBox.xMax = cHigh;

  // set the bottom coordinate
  cLow = current.yMin;
  cHigh = current.yMax;
  nLow = next.yMin;
  nHigh = next.yMax;
  HBBoxShortToLeft(cLow, cHigh, nLow, nHigh);
  newBox.yMin = cLow;
  newBox.yMax = cHigh;
  return newBox;
}

BBox HBBoxShrinkX2Center(const BBox& current, const double& scale) {
  double width, left, newWidth;
  width = current.xMax - current.xMin;
  newWidth = width * scale;
  left = current.xMin + ((width - newWidth) / 2);
  return BBox(left, current.yMin, left + newWidth, current.yMax);
}

BBox HBBoxShrinkY2Center(const BBox& current, const double& scale) {
  double height, top, newHeight;
  height = current.yMax - current.yMin;
  newHeight = height * scale;
  top = current.yMax - ((height - newHeight) / 2);
  return BBox(current.xMin, top - newHeight, current.xMax, top);
}

BBox HBBoxShrink2Right(const BBox& current, const double& scale) {
  double width, left;
  width = (current.xMax - current.xMin) * scale;
  left = current.xMax - width;
  return BBox(left, current.yMin, current.xMax, current.yMax);
}

BBox HBBoxShrink2Bottom(const BBox& current, const double& scale) {
  double height, top;
  height = (current.yMax - current.yMin) * scale;
  top = current.yMin + height;
  return BBox(current.xMin, current.yMin, current.xMax, top);
}

BBox HBBoxShrink2Left(const BBox& current, const double& scale) {
  double width, right;
  width = (current.xMax - current.xMin) * scale;
  right = current.xMin + width;
  return BBox(current.xMin, current.yMin, right, current.yMax);
}

BBox HBBoxShrink2Top(const BBox& current, const double& scale) {
  double height, bottom;
  height = (current.yMax - current.yMin) * scale;
  bottom = current.yMax - height;
  return BBox(current.xMin, bottom, current.xMax, current.yMax);
}
