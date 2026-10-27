// Regression for clipping sparse gain buckets, including an empty negative side.
#include <iostream>
#include <vector>

#include "gainBucket.h"
#include "gainBucketHL.h"

static void checkClip(const std::vector<int>& gains, unsigned maxGain) {
  BucketArray dense;
  BucketArrayHL sparse;
  dense.setMaxGain(maxGain);
  sparse.setMaxGain(maxGain);
  std::vector<SVGainElement> denseElements(gains.size());
  std::vector<SVGainElement> sparseElements(gains.size());
  for (unsigned i = 0; i < gains.size(); ++i) {
    denseElements[i]._gain = sparseElements[i]._gain = gains[i];
    dense.push(denseElements[i]);
    sparse.push(sparseElements[i]);
  }

  dense.setupForClip();
  sparse.setupForClip();
  SVGainElement* expected = dense[maxGain]._next;
  SVGainElement* actual = sparse[maxGain]._next;
  for (unsigned i = 0; i < gains.size(); ++i) {
    abkfatal(expected && actual, "Clipping lost a gain element");
    const unsigned index = expected - denseElements.data();
    abkfatal(actual == &sparseElements[index], "Clipping changed move order");
    abkfatal(actual->getGain() == 0 && actual->getRealGain() == gains[index],
             "Clipping did not preserve the real gain");
    expected = expected->_next;
    actual = actual->_next;
  }
  abkfatal(!expected && !actual, "Clipping duplicated an element or made a cycle");
  abkfatal(sparse.isConsistent(), "Clipping broke gain-list links");
  for (unsigned i = 0; i < gains.size(); ++i)
    if (gains[i] != 0)
      abkfatal(sparse[gains[i] + maxGain].isEmpty(),
               "Clipping left an element in its old bucket");
}

int main() {
  checkClip({}, 10);
  checkClip({0}, 10);
  checkClip({0, 0}, 0);
  checkClip({1, 3, 1}, 10);
  checkClip({0, 3, 0, 1}, 10);
  checkClip({-1, -3, -1}, 10);
  checkClip({-3, 0, 2, -1, 2}, 10);
  std::cout << "Sparse gain clipping checks passed" << std::endl;
}
