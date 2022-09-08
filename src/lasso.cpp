#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <math.h>
#include <H5Cpp.h>

#define PI 3.14159265
/*
const std::string FileName("results.h5");
const std::string DatasetName("PersonalInformation");
const std::string member_rho("rho ");
const std::string member_u("u");
const std::string member_v("v");

*/

const int nDim = 2; // TODO: EXAMPLE!
const int nDist = 9; // TODO: EXAMPLE!
const float w [9] = {4.0 / 9.0, 
                     1.0 / 9.0,
                     1.0 / 9.0,
                     1.0 / 9.0,
                     1.0 / 9.0,
                     1.0 / 36.0,
                     1.0 / 36.0,
                     1.0 / 36.0,
                     1.0 / 36.0};

const float c [2][9] = { {0. , 1. , 0. , -1. , 0. , 1. , -1. , -1. , 1.},
                         {0. , 0. , 1. , 0. , -1. , 1. , 1. , -1. , -1.} };

const float cs = sqrt(1. / 3.);
const float csSq = 1. / 3.;

const float nu = 0.5;
const float tau_s = (6. * nu + 1.) / 2.;

namespace parameters {
  const float p_left = 0.1;
  const float p_right = 0.8;
  const bool sourceTerms = true;
  struct monopole {
    const float amplitude = 0.01;
    const float frequency = 0.25;
  } monopole1;
}

struct cell {
  float rho;
  float u [2]; // TODO: EXAMPLE!
  float dist [9]; // TODO: EXAMPLE!
  float eqDist [9]; // TODO: EXAMPLE!
  float fstar [9]; // TODO: EXAMPLE!;
  bool isBoundary = false;
  float x;
  float y;
};

void computeEquilibriumDistribution(const float rho, const float* u, float * eqDist) {
  for (int i = 0; i < nDist; i++) {
    eqDist[i] = w[i] * rho * ( 1. +  (u[0] * c [0][i] + u[1] * c[1][i]) / csSq +  (u[0] * c [0][i] + u[1] * c[1][i]) * (u[0] * c [0][i] + u[1] * c[1][i])  / (2. * csSq *  csSq) 
                                                                                - (u[0] * u[0] + u[1] * u[1]) / (2. * csSq) );
    // std::cout << i << "\t\t"<< eqDist[i] << std::endl;
  }
}

int main(int argc, char** argv) {
  // Welcome message to relieve the stress of the users
  std::cout << "Welcome to LASSO!" << std::endl;
  std::cout << "Your 2D flow solver using the lattice Boltzmann method." << std::endl;

  std::cout << "Created by Benyamin Krisna" << std::endl;

  const int N = 32; // TODO: EXAMPLE!


  // Initialization
  cell Grid [N][N];
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      /*
      if (i * j == 0 || j == N-1 || i == N-1 ) {
        Grid[i][j].isBoundary = true;
      }
      */
      if (Grid[i][j].isBoundary == false) {
        Grid[i][j].rho = 1.0;
        Grid[i][j].u[0] = 0.0;
        Grid[i][j].u[1] = 0.0;
        // std::cout << "Grid location: " << i << "\t\t" << j << std::endl;
        computeEquilibriumDistribution(Grid[i][j].rho, Grid[i][j].u, Grid[i][j].eqDist);
        for (int dist_count = 0; dist_count < nDist; dist_count++) {
          Grid[i][j].dist[dist_count] = Grid[i][j].eqDist[dist_count]; // set the distribution to equilibrium distribution
        }
      } else { // is a boundary cell

      }
    }
  }

  // Initialize boundary cells
  // Top wall
  for (int i = 0; i < N; i++) {
    Grid[i][N-1].u[0] = 0.;
    Grid[i][N-1].u[1] = 0.;
  }

  // Bottom wall
  for (int i = 0; i < N; i++) {
    Grid[i][0].u[0] = 0.;
    Grid[i][0].u[1] = 0.;
  }

  // Time step algorithm
  // 1. Compute macroscopic moments from the distributions
  // 2. Obtain equilibrium distributions
  // 3. Write macroscopic  fields
  // 4. Perform collision (relaxation)
  // 5. Perform streaming (propagation)
  // 6. Increase the time step
  // 7. Repeat!

  int globalTimeStep = 0;
  int tMax = 2000;

  while (globalTimeStep < tMax) {
    // 1. Compute macroscopic moments from the distribution
    for (int i = 0; i < N; i++) {
      for (int j = 0; j < N; j++) {
        if (!Grid[i][j].isBoundary) {
          float a_rho = 0.;
          float a_u = 0.;
          float a_v = 0.;
          for (int dist_count = 0; dist_count < nDist; dist_count++) {
            a_rho += Grid[i][j].dist[dist_count];
            a_u += ( c[0][dist_count] * Grid[i][j].dist[dist_count]) ;
            a_v += ( c[1][dist_count] * Grid[i][j].dist[dist_count]) ;
          }
          a_u /= a_rho;
          a_v /= a_rho;
          Grid[i][j].rho = a_rho;
          Grid[i][j].u[0] = a_u;
          Grid[i][j].u[1] = a_v;
          // std::cout << i << "\t\t" << j << "\t\t" << a_rho << std::endl;
        }
      }
    }

    // 2. Obtain equilibrium distributions
    for (int i = 0; i < N; i++) {
      for (int j = 0; j < N; j++) {
        computeEquilibriumDistribution(Grid[i][j].rho, Grid[i][j].u, Grid[i][j].eqDist);
      }
    }

    // 3. Write macroscopic fields

    // 4. Perform collision  (relaxation)
    for (int i = 0; i < N; i++) {
      for (int j = 0; j < N; j++) {
        for (int dist_count = 0; dist_count < nDist; dist_count++) {
          Grid[i][j].fstar[dist_count] = Grid[i][j].dist[dist_count] * ( 1. - 1. / tau_s ) + Grid[i][j].eqDist[dist_count] / tau_s;
        }
      }
    }
    
    // 4b. Any source terms / forces?
    if (parameters::sourceTerms) {
      //TODO: very  hard coded!
      for (int dist_count = 0; dist_count < nDist; dist_count++) {
        Grid[N/2][N/2].fstar[dist_count] += w[dist_count] * parameters::monopole1.amplitude * sin(2.0 * PI * parameters::monopole1.frequency);
      }
    }

    // 5. Propagation step
    for (int i = 1; i < N-1; i++) {
      for (int j = 1; j < N-1; j++) {
        for (int dist_count = 0; dist_count < nDist; dist_count++) {
          Grid[i + (int)c[0][dist_count]][j + (int)c[1][dist_count]].dist[dist_count] = Grid[i][j].fstar[dist_count];
        }
      }
    }

    /*
    // 5b. Propagation for boundary cells
    for (int i = 1; i < N-1; i++) {
      // Bottom wall
      Grid[i][0].dist[2] = Grid[i][0].fstar[4];
      Grid[i][0].dist[5] = Grid[i][0].fstar[7];
      Grid[i][0].dist[6] = Grid[i][0].fstar[8];

      // Top wall
      Grid[i][N-1].dist[4] = Grid[i][N-1].fstar[2];
      Grid[i][N-1].dist[7] = Grid[i][N-1].fstar[5];
      Grid[i][N-1].dist[8] = Grid[i][N-1].fstar[6];

        // inlet (left boundary)
      const float rho_left = parameters::p_left / csSq; 
      const float rho_righ = parameters::p_right / csSq;
      for (int dist_count = 0; dist_count < nDist; dist_count++) {
        //TODO: phrase?
        //TODO: define u_w and rho_w
        Grid[0][i].dist[1] = - Grid[0][i].fstar[3] + 2 * w[dist_count] * rho_w * ( 1. + (u[0] * c [0][i] + u[1] * c[1][i]) * (u[0] * c [0][i] + u[1] * c[1][i])  / (2. * csSq *  csSq) - phrase );
      }
    }
    */

    // 6. Increase the time step
    globalTimeStep++;
    
  }

  // TODO: hard-coded output replacement?
  std::ofstream output_file;
  output_file.open("results.dat");
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      output_file << i <<  '\t' << j << '\t' << (Grid[i][j].rho - 1.0) << std::endl;
    }
    output_file << std::endl;
  }

  std::ofstream sample_file;
  samplsample.open("sampling.dat");
  for (int i = 0; i < N; i++) {
    sample_file << i << '\t' << (Grid[i][N/2].rho - 1.0) << std::endl;
  }

}
