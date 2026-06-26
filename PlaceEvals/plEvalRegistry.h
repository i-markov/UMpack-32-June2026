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

// created by Stefanus Mantik on 06/19/98

#include <iostream>

#include <string>
#include <vector>

class HGraphFixed;
class Placement;
class LayoutBBoxes;
class PlaceEvalXFace;

#ifndef _PLACEEVALREG_H_
#define _PLACEEVALREG_H_

class PlaceEvalType {
 public:
  enum Type {
    BBoxPlEval,
    BBoxPlEvalWt,
    BBoxWCheng,
    BBoxWRSMT,
    HBBox,
    HBBoxWCheng,
    HBBoxWRSMT,
    HBBox0,
    HBBox0wCheng,
    HBBox0wRSMT,
    HBBoxN,
    HBBoxNwCheng,
    HBBoxNwRSMT,
    HBBox6,
    HBBox6wRSMT,
    OneSteiner
  };
  // single-cell move stucture
  // is the default for names,
  // so is "with net tallies"
 protected:
  Type _type;
  static char strbuf[];

 public:
  PlaceEvalType(Type type = BBoxPlEval) : _type(type){};
  PlaceEvalType(int argc, const char* argv[]);  // catches -placeEval <name>

  operator unsigned() const { return _type; }

  double computeCost(const HGraphFixed&, const Placement&);
  double computeCost(const HGraphFixed&, const LayoutBBoxes&,
                     std::vector<unsigned>&);
  double computeCost(const HGraphFixed&, const LayoutBBoxes&, const Placement&);

  void costPerNet(const HGraphFixed&, const Placement&);
  void costPerNet(const HGraphFixed&, const LayoutBBoxes&,
                  std::vector<unsigned>&);
  void costPerNet(const HGraphFixed&, const LayoutBBoxes&, const Placement&);

  void comparison(const HGraphFixed&, const LayoutBBoxes&, const Placement&,
                  const char* outName);

  PlaceEvalXFace* init(const HGraphFixed&, const LayoutBBoxes&,
                       const Placement&);

  operator const char*() const;
};

std::ostream& operator<<(std::ostream& os, const PlaceEvalType& pet);

#endif
