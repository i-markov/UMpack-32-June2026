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
template <typename ResourceT, typename ResourceDiffT, typename ResourceZeroT>
Generic2DResourceMap<ResourceT, ResourceDiffT,
                     ResourceZeroT>::Generic2DResourceMap() {}

template <typename ResourceT, typename ResourceDiffT, typename ResourceZeroT>
Generic2DResourceMap<ResourceT, ResourceDiffT,
                     ResourceZeroT>::~Generic2DResourceMap() {}

// input functions

// will overwrite current map with a new default constructed one of the
// specified size
template <typename ResourceT, typename ResourceDiffT, typename ResourceZeroT>
void Generic2DResourceMap<ResourceT, ResourceDiffT, ResourceZeroT>::setNumTiles(
    unsigned numHorizTiles, unsigned numVertTiles) {
  _numHorizTiles = numHorizTiles;
  _numVertTiles = numVertTiles;

  // not sure how to do this part just using my typedefs, just going to use
  // vector for now, sorry.
  _demandMap = std::vector<std::vector<ResourceT> >(
      numHorizTiles, std::vector<ResourceT>(numVertTiles, ResourceZeroT()()));
  _supplyMap = std::vector<std::vector<ResourceT> >(
      numHorizTiles, std::vector<ResourceT>(numVertTiles, ResourceZeroT()()));
}

template <typename ResourceT, typename ResourceDiffT, typename ResourceZeroT>
void Generic2DResourceMap<ResourceT, ResourceDiffT, ResourceZeroT>::setTileSize(
    double tileWidth, double tileHeight) {
  _tileWidth = tileWidth;
  _tileHeight = tileHeight;
}

template <typename ResourceT, typename ResourceDiffT, typename ResourceZeroT>
void Generic2DResourceMap<ResourceT, ResourceDiffT, ResourceZeroT>::
    setLayoutBBox(const BBox& layoutBBox) {
  _layoutBBox = layoutBBox;
}

template <typename ResourceT, typename ResourceDiffT, typename ResourceZeroT>
void Generic2DResourceMap<ResourceT, ResourceDiffT, ResourceZeroT>::
    setDemandMap(const ResourceMap& demandMap) {
  _demandMap = demandMap;
}

template <typename ResourceT, typename ResourceDiffT, typename ResourceZeroT>
void Generic2DResourceMap<ResourceT, ResourceDiffT, ResourceZeroT>::
    setDemandValue(unsigned x, unsigned y, const ResourceT& value) {
  abkassert(x < _demandMap.size(), "Range error in Generic2DResourceMap");
  abkassert(0 <= x, "Range error in Generic2DResourceMap");
  abkassert(y < _demandMap[x].size(), "Range error in Generic2DResourceMap");
  abkassert(0 <= y, "Range error in Generic2DResourceMap");

  _demandMap[x][y] = value;
}

template <typename ResourceT, typename ResourceDiffT, typename ResourceZeroT>
void Generic2DResourceMap<ResourceT, ResourceDiffT, ResourceZeroT>::
    setSupplyMap(const ResourceMap& supplyMap) {
  _supplyMap = supplyMap;
}

template <typename ResourceT, typename ResourceDiffT, typename ResourceZeroT>
void Generic2DResourceMap<ResourceT, ResourceDiffT, ResourceZeroT>::
    setSupplyValue(unsigned x, unsigned y, const ResourceT& value) {
  abkassert(x < _supplyMap.size(), "Range error in Generic2DResourceMap");
  abkassert(0 <= x, "Range error in Generic2DResourceMap");
  abkassert(y < _supplyMap[x].size(), "Range error in Generic2DResourceMap");
  abkassert(0 <= y, "Range error in Generic2DResourceMap");

  _supplyMap[x][y] = value;
}

template <typename ResourceT, typename ResourceDiffT, typename ResourceZeroT>
BBox Generic2DResourceMap<ResourceT, ResourceDiffT, ResourceZeroT>::getTileBBox(
    unsigned i, unsigned j) const {
  double left = i * _tileWidth + _layoutBBox.xMin;
  double bottom = j * _tileHeight + _layoutBBox.yMin;
  double right = left + _tileWidth;
  double top = bottom + _tileHeight;

  BBox rval(left, bottom, right, top);
  return rval;
}

// inspectors
template <typename ResourceT, typename ResourceDiffT, typename ResourceZeroT>
ResourceT
Generic2DResourceMap<ResourceT, ResourceDiffT, ResourceZeroT>::getDemand(
    const BBox& bbox) const {
  // TODO
  abkfatal(false, "Function unimplemented");
  return ResourceZeroT()();
}

template <typename ResourceT, typename ResourceDiffT, typename ResourceZeroT>
ResourceT Generic2DResourceMap<ResourceT, ResourceDiffT,
                               ResourceZeroT>::getDemand(unsigned x,
                                                         unsigned y) const {
  abkassert(x < _demandMap.size(), "Range error in Generic2DResourceMap");
  abkassert(0 <= x, "Range error in Generic2DResourceMap");
  abkassert(y < _demandMap[x].size(), "Range error in Generic2DResourceMap");
  abkassert(0 <= y, "Range error in Generic2DResourceMap");

  return _demandMap[x][y];
}

template <typename ResourceT, typename ResourceDiffT, typename ResourceZeroT>
ResourceT
Generic2DResourceMap<ResourceT, ResourceDiffT, ResourceZeroT>::getSupply(
    const BBox& bbox) const {
  // TODO
  abkfatal(false, "Function unimplemented");
  return ResourceZeroT()();
}

template <typename ResourceT, typename ResourceDiffT, typename ResourceZeroT>
ResourceT Generic2DResourceMap<ResourceT, ResourceDiffT,
                               ResourceZeroT>::getSupply(unsigned x,
                                                         unsigned y) const {
  abkassert(x < _supplyMap.size(), "Range error in Generic2DResourceMap");
  abkassert(0 <= x, "Range error in Generic2DResourceMap");
  abkassert(y < _supplyMap[x].size(), "Range error in Generic2DResourceMap");
  abkassert(0 <= y, "Range error in Generic2DResourceMap");

  return _supplyMap[x][y];
}

template <typename ResourceT, typename ResourceDiffT, typename ResourceZeroT>
double Generic2DResourceMap<ResourceT, ResourceDiffT,
                            ResourceZeroT>::getTotalOverfill(void) const {
  double overfill = 0.;

  for (unsigned i = 0; i < _numHorizTiles; ++i) {
    for (unsigned j = 0; j < _numVertTiles; ++j) {
      overfill +=
          std::max(0., ResourceDiffT()(getDemand(i, j), getSupply(i, j)));
    }
  }

  return overfill;
}

template <typename ResourceT, typename ResourceDiffT, typename ResourceZeroT>
double Generic2DResourceMap<ResourceT, ResourceDiffT,
                            ResourceZeroT>::getMaxOverfill(void) const {
  double maxOverfill = 0.;

  for (unsigned i = 0; i < _numHorizTiles; ++i) {
    for (unsigned j = 0; j < _numVertTiles; ++j) {
      maxOverfill = std::max(
          maxOverfill,
          std::max(0., ResourceDiffT()(getDemand(i, j), getSupply(i, j))));
    }
  }

  return maxOverfill;
}

template <typename ResourceT, typename ResourceDiffT, typename ResourceZeroT>
unsigned Generic2DResourceMap<ResourceT, ResourceDiffT,
                              ResourceZeroT>::getNumOverfullTiles(void) const {
  unsigned num = 0;

  for (unsigned i = 0; i < _numHorizTiles; ++i) {
    for (unsigned j = 0; j < _numVertTiles; ++j) {
      if (ResourceDiffT()(getDemand(i, j), getSupply(i, j)) < 0) {
        ++num;
      }
    }
  }

  return num;
}

template <typename ResourceT, typename ResourceDiffT, typename ResourceZeroT>
void Generic2DResourceMap<ResourceT, ResourceDiffT,
                          ResourceZeroT>::printDemandMap(void) const {
  for (unsigned rj = 1; rj <= _numVertTiles; ++rj) {
    unsigned j = _numVertTiles - rj;
    for (unsigned i = 0; i < _numHorizTiles; ++i) {
      std::cout << getDemand(i, j);
    }
    std::cout << std::endl;
  }
}

template <typename ResourceT, typename ResourceDiffT, typename ResourceZeroT>
void Generic2DResourceMap<ResourceT, ResourceDiffT,
                          ResourceZeroT>::printSupplyMap(void) const {
  for (unsigned rj = 1; rj <= _numVertTiles; ++rj) {
    unsigned j = _numVertTiles - rj;
    for (unsigned i = 0; i < _numHorizTiles; ++i) {
      std::cout << getSupply(i, j);
    }
    std::cout << std::endl;
  }
}

// TODO, fix this function to be generic and work with ResourceT
// Ploting functions in various formats
template <typename ResourceT, typename ResourceDiffT, typename ResourceZeroT>
void Generic2DResourceMap<ResourceT, ResourceDiffT, ResourceZeroT>::plotGP(
    const std::string& baseFileName) const {
  std::string fname = baseFileName + ".plt";

  std::cout << "Saving " << fname << std::endl;

  std::ofstream gpFile(fname.c_str());

  gpFile << "#Use this file as a script for gnuplot" << std::endl
         << "#(See http://www.gnuplot.info/ for details)" << std::endl
         << "#" << TimeStamp() << User() << Platform() << std::endl
         << std::endl
         << "set title ' " << baseFileName << std::endl
         << std::endl
         << "set nokey" << std::endl
         << "#   Uncomment these two lines starting with \"set\"" << std::endl
         << "#   to save an EPS file for inclusion into a latex document"
         << std::endl
         << "# set terminal postscript eps color solid 10" << std::endl
         << "# set output \"" << baseFileName << ".eps\"" << std::endl
         << std::endl
         << "#   Uncomment these two lines starting with \"set\"" << std::endl
         << "#   to save a PS file for printing" << std::endl
         << "# set terminal postscript portrait color solid 8" << std::endl
         << "# set output \"" << baseFileName << ".ps\"" << std::endl
         << std::endl
         << std::endl
         << "set view 55,45 " << std::endl
         << std::endl
         << std::endl
         << "splot '-' w l " << std::endl
         << std::endl;

  for (unsigned i = 0; i < _numHorizTiles; ++i) {
    for (unsigned j = 0; j < _numVertTiles; ++j) {
      double regXMin = _layoutBBox.xMin + static_cast<double>(i) * _tileWidth;
      double regXMax =
          _layoutBBox.xMin + static_cast<double>(i + 1) * _tileWidth;
      double regYMin = _layoutBBox.yMin + static_cast<double>(j) * _tileHeight;
      double regYMax =
          _layoutBBox.yMin + static_cast<double>(j + 1) * _tileHeight;

      gpFile << regXMin << " " << regYMin << " "
             << getMagnitude(getDemand(i, j)) << std::endl;
      gpFile << regXMax << " " << regYMin << " "
             << getMagnitude(getDemand(i, j)) << std::endl;
      gpFile << regXMax << " " << regYMax << " "
             << getMagnitude(getDemand(i, j)) << std::endl;
      gpFile << regXMin << " " << regYMax << " "
             << getMagnitude(getDemand(i, j)) << std::endl;
      gpFile << regXMin << " " << regYMin << " "
             << getMagnitude(getDemand(i, j)) << std::endl;
    }
  }

  gpFile << "EOF" << std::endl;
  gpFile << "pause -1 'Press any key' " << std::endl;
  gpFile.close();
}

// TODO, fix this function to be generic and work with ResourceT
template <typename ResourceT, typename ResourceDiffT, typename ResourceZeroT>
void Generic2DResourceMap<ResourceT, ResourceDiffT, ResourceZeroT>::plotMatlab(
    const std::string& baseFileName) const {
  std::string fname = baseFileName + ".m";

  std::cout << "Saving " << fname << std::endl;

  std::ofstream gpFile(fname.c_str());

  gpFile << "X = [ " << std::endl;
  for (unsigned i = 0; i < _numHorizTiles; ++i) {
    double regXMin = _layoutBBox.xMin + static_cast<double>(i) * _tileWidth;
    double regXMax = _layoutBBox.xMin + static_cast<double>(i + 1) * _tileWidth;
    gpFile << 0.5 * (regXMin + regXMax) << " ; " << std::endl;
  }
  gpFile << "] ;" << std::endl << std::endl;

  gpFile << "Y = [ " << std::endl;
  for (unsigned i = 0; i < _numVertTiles; ++i) {
    double regYMin = _layoutBBox.yMin + static_cast<double>(i) * _tileHeight;
    double regYMax = _layoutBBox.xMin + static_cast<double>(i + 1) * _tileWidth;
    gpFile << 0.5 * (regYMin + regYMax) << " ; " << std::endl;
  }
  gpFile << "] ;" << std::endl << std::endl;

  gpFile << "Z = [" << std::endl << std::endl;
  for (unsigned j = 0; j < _numVertTiles; ++j) {
    for (unsigned i = 0; i < _numHorizTiles; ++i) {
      gpFile << getMagnitude(getDemand(i, j)) << " ";
    }
    gpFile << " ; " << std::endl;
  }
  gpFile << "] ;" << std::endl << std::endl;

  gpFile << "surf(X,Y,Z)" << std::endl << std::endl;
  gpFile << "colormap hot" << std::endl;

  gpFile << "xlabel('X Axis')" << std::endl;
  gpFile << "ylabel('Y Axis')" << std::endl;
  gpFile << "zlabel('Congestion')" << std::endl;
  gpFile << "title('Total Congestion')" << std::endl;

  gpFile << "view(0,90)" << std::endl;

  gpFile.close();
}

// TODO, fix this function to be generic and work with ResourceT
template <typename ResourceT, typename ResourceDiffT, typename ResourceZeroT>
void Generic2DResourceMap<ResourceT, ResourceDiffT, ResourceZeroT>::plotXPM(
    const std::string& baseFileName, const ResourceT& maxResourceVal) const {
  unsigned numColors = 64;
  Palette palette = getPalette();
  std::vector<std::string> ColorMap = palette.first;
  std::vector<std::string> Colors = palette.second;

  std::string fname = baseFileName + ".xpm";

  std::cout << "Saving " << fname << std::endl;

  std::ofstream xpmFile(fname.c_str());

  ResourceT maxResourceValue = maxResourceVal;

  unsigned xBlowFactor =
      unsigned(ceil(400.0 / static_cast<double>(_numHorizTiles)));
  unsigned yBlowFactor =
      unsigned(ceil(400.0 / static_cast<double>(_numVertTiles)));

  xpmFile << "/* XPM */" << std::endl;
  xpmFile << "static char *congestion[] = {" << std::endl;
  xpmFile << "/* columns rows colors chars-per-pixel */" << std::endl;
  xpmFile << "\"" << _numHorizTiles * xBlowFactor << " "
          << _numVertTiles * yBlowFactor << " " << numColors << " 1\","
          << std::endl;
  for (unsigned i = 0; i < numColors; ++i)
    xpmFile << "\"" << ColorMap[i] << "\"," << std::endl;
  xpmFile << "/* pixels */" << std::endl;

  std::vector<std::vector<double> > image;
  std::vector<std::pair<unsigned, unsigned> > maxLocs;
  for (unsigned rj = 1; rj <= _numVertTiles; ++rj) {
    unsigned j = _numVertTiles - rj;
    std::vector<double> horizLine;
    for (unsigned i = 0; i < _numHorizTiles; ++i) {
      ResourceT consumption = getDemand(i, j);
      if (ResourceDiffT()(consumption, maxResourceValue) > 0.) {
        maxLocs.clear();
        maxLocs.push_back(std::make_pair(i, rj - 1));
        maxResourceValue = consumption;
      } else if (consumption == maxResourceValue) {
        maxLocs.push_back(std::make_pair(i, rj - 1));
      }
      horizLine.push_back(getMagnitude(consumption));
    }
    image.push_back(horizLine);
  }

  bool drawBoxAroundHighest = false;
  if (drawBoxAroundHighest) {
    for (unsigned radius = 30; radius <= 34; ++radius) {
      for (unsigned i = 0; i < maxLocs.size(); ++i) {
        unsigned x = maxLocs[i].first;
        unsigned y = maxLocs[i].second;

        // draw left side
        if (x >= radius) {
          unsigned j = y;
          for (unsigned k = 0; k <= radius; ++k) {
            image[j][x - radius] = getMagnitude(maxResourceValue);
            if (j == 0) break;
            --j;
          }
          j = y;
          for (unsigned k = 0; k <= radius; ++k) {
            image[j][x - radius] = getMagnitude(maxResourceValue);
            ++j;
            if (j == _numVertTiles) break;
          }
        }

        // draw right side
        if (_numHorizTiles - x > radius) {
          unsigned j = y;
          for (unsigned k = 0; k <= radius; ++k) {
            image[j][x + radius] = getMagnitude(maxResourceValue);
            if (j == 0) break;
            --j;
          }
          j = y;
          for (unsigned k = 0; k <= radius; ++k) {
            image[j][x + radius] = getMagnitude(maxResourceValue);
            ++j;
            if (j == _numVertTiles) break;
          }
        }

        // draw bottom
        if (y >= radius) {
          unsigned j = x;
          for (unsigned k = 0; k <= radius; ++k) {
            image[y - radius][j] = getMagnitude(maxResourceValue);
            if (j == 0) break;
            --j;
          }
          j = x;
          for (unsigned k = 0; k <= radius; ++k) {
            image[y - radius][j] = getMagnitude(maxResourceValue);
            ++j;
            if (j == _numHorizTiles) break;
          }
        }

        // draw top
        if (_numVertTiles - y > radius) {
          unsigned j = x;
          for (unsigned k = 0; k <= radius; ++k) {
            image[y + radius][j] = getMagnitude(maxResourceValue);
            if (j == 0) break;
            --j;
          }
          j = x;
          for (unsigned k = 0; k <= radius; ++k) {
            image[y + radius][j] = getMagnitude(maxResourceValue);
            ++j;
            if (j == _numHorizTiles) break;
          }
        }
      }
    }
  }

  std::cout << "The maximum resource value is " << maxResourceValue
            << std::endl;

  for (unsigned j = 0; j < _numVertTiles; ++j) {
    for (unsigned k = 0; k < yBlowFactor; ++k) {
      xpmFile << "\"";
      for (unsigned i = 0; i < _numHorizTiles; ++i) {
        double conversion = (numColors - 1) / getMagnitude(maxResourceValue);
        unsigned index = 0;
        if (equalDouble(image[j][i],
                        0.0))  // to specialize the case where demand==supply==0
          index = 0;
        else
          index = unsigned(floor(image[j][i] * conversion));

        for (unsigned l = 0; l < xBlowFactor; ++l) {
          xpmFile << Colors[index];
        }
      }
      if (j == _numVertTiles - 1 && k == yBlowFactor - 1) {
        xpmFile << "\"" << std::endl;
      } else {
        xpmFile << "\"," << std::endl;
      }
    }
  }
  xpmFile << "};" << std::endl;

  xpmFile.close();
}

// get the color palette to use in plotting
template <typename ResourceT, typename ResourceDiffT, typename ResourceZeroT>
Palette Generic2DResourceMap<ResourceT, ResourceDiffT,
                             ResourceZeroT>::getPalette(void) const {
  return PaletteBox::Hot();
}

template <typename ResourceT, typename ResourceDiffT, typename ResourceZeroT>
double
Generic2DResourceMap<ResourceT, ResourceDiffT, ResourceZeroT>::getMagnitude(
    const ResourceT& r) const {
  return ResourceDiffT()(r, ResourceZeroT()());
}

template <typename ResourceT, typename ResourceDiffT, typename ResourceZeroT>
double Generic2DResourceMap<ResourceT, ResourceDiffT, ResourceZeroT>::
    estimateOverfullnessChangeFromCellMove(unsigned cellId,
                                           const Point& newPos) const {
  return 0.;
}
