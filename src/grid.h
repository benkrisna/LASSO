#ifndef __GRID__
#define __GRID__

#include <array>

/*
typedef struct {
  float rho;
  // float u [2]; // TODO: EXAMPLE!
  // float dist [9]; // TODO: EXAMPLE!
  // float eqDist [9]; // TODO: EXAMPLE!
  // float fstar [9];
  float* u;
  float* dist;
  float* eqDist;
  float* fstar;
  float x;
  float y;

  bool isBoundary = false;
} cell;
*/

class cell {
  public:
    float rho;
    // float u [2]; // TODO: EXAMPLE!
    // float dist [9]; // TODO: EXAMPLE!
    // float eqDist [9]; // TODO: EXAMPLE!
    // float fstar [9];
    float* u;
    float* dist;
    float* eqDist;
    float* fstar;
    float x;
    float y;

    bool isBoundary = false;
    float& Rho() {return rho;};
    float& X() {return x;};
    float& Y() {return y;};
};

class LassoGrid {
  public:
    LassoGrid(const int m_N, const int m_nDim, const int m_nDist);
    int gridPoints() const;
    cell& getCell(int* index);
    void initializeGrid();

  private:
    const int N;
    const int nDim;
    const int nDist;
    cell* Grid;
};

#endif
