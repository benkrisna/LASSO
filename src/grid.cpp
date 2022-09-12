#include "grid.h"
#include "math.h"

#include <iostream>
#include <array>

LassoGrid::LassoGrid(const int m_N, const int m_nDim, const int m_nDist) : N {m_N}, nDim {m_nDim}, nDist {m_nDist} {
  // Allocate Grid
  const int numberCells = pow(N, nDim);
  Grid = new cell[numberCells];
  std::cout << numberCells << std::endl;
  for (int a = 0; a < numberCells; a++) {
    Grid[a].u = new float [nDim];
    Grid[a].dist = new float [nDist];
    Grid[a].eqDist = new float [nDist];
    Grid[a].fstar = new float [nDist];
  }
}

int LassoGrid::gridPoints() const {
  return N;
}

cell& LassoGrid::getCell(int* index) {
  int thisIndex = 0;
  for (int d = nDim - 1; d >= 0; d--) {
    thisIndex += pow(N, d) * index[d];
  }
  return Grid[thisIndex];
}
