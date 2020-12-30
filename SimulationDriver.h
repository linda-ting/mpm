#include <Eigen/Sparse>
#include <unsupported/Eigen/IterativeSolvers>

#include <sys/stat.h>
#include <iostream>
#include "MaterialPointSystem.h"
#include <functional>

class SimulationDriver{
public:
    MaterialPointSystem ms = MaterialPointSystem(100000);
    float dt;
    Eigen::Matrix<float,3,1> gravity;
    std::string test;

    SimulationDriver()
    : dt((float)0.00001)
    {
      gravity.setZero();
      gravity(1) = -9.8;
    }

    void run(const int max_frame)
    {
      float accumulate_t = 0;
      mkdir("output/", 0777);
      std::string output_folder = "output/" + test + "_" + std::to_string(ms.youngs_modulus);
      mkdir(output_folder.c_str(), 0777);
      std::string filename = output_folder + "/" + std::to_string(0) + ".poly";
      ms.dumpPoly(filename);
      for(int frame=1; frame<=max_frame; frame++) {
          std::cout << "Frame " << frame << std::endl;
          int N_substeps = (int)(((float)1/24)/dt);
          for (int step = 1; step <= N_substeps; step++) {
              std::cout << "Step " << step << std::endl;
              advanceOneStep();
              accumulate_t += dt;
          }
          mkdir("output/", 0777);
          std::string output_folder = "output/" + test + "_" + std::to_string(ms.youngs_modulus);
          mkdir(output_folder.c_str(), 0777);
          std::string filename = output_folder + "/" + std::to_string(frame) + ".poly";
          ms.dumpPoly(filename);
          std::cout << std::endl;
      }
    }

    void advanceOneStep()
    {
      int x = ms.grid->dims[0];
      int y = ms.grid->dims[1];
      int z = ms.grid->dims[2];
      float h = ms.grid->h;      
      int n = ms.n;

      // clear data; zero mass, velocity, force
      for (int i = 0; i < x; i++) {
        for (int j = 0; j < y; j++) {
          for (int k = 0; k < z; k++) {
            ms.grid->m[i][j][k] = 0.f;
            ms.grid->v[i][j][k] = Eigen::Matrix<float,3,1>::Zero();
            ms.grid->f[i][j][k] = Eigen::Matrix<float,3,1>::Zero();
          }
        }
      }
      std::cout << "cleared data" << std::endl;

      // transfer mass and momentum from particles to grid
      for (int m = 0; m < n; m++) {
        Particle* p = ms.particles[m];
        Eigen::Matrix<float,3,1> p_pos = p->x;

        int xb = int((p_pos[0] - 0.5f * h) / h);
        int yb = int((p_pos[1] - 0.5f * h) / h);
        int zb = int((p_pos[2] - 0.5f * h) / h);
        xb = std::max(0, xb);
        yb = std::max(0, yb);
        zb = std::max(0, zb);
        
        for (int i = xb; i <= xb + 2 && i < x; i++) {
          for (int j = yb; j <= yb + 2 && j < y; j++) {
            for (int k = zb; k <= zb + 2 && k < z; k++) {
              Eigen::Matrix<float,3,1> n_pos = ms.grid->getPosition(i, j, k);
              float w = ms.omega(p_pos, n_pos);
              ms.grid->m[i][j][k] += p->m * w;
              ms.grid->v[i][j][k] += p->m * p->v * w;
            }
          }
        }
      }
      std::cout << "transfered from particles to grid" << std::endl;

      // calculate each grid node's velocity
      for (int i = 0; i < x; i++) {
        for (int j = 0; j < y; j++) {
          for (int k = 0; k < z; k++) {
            if (ms.grid->m[i][j][k] > 0.f) {
              ms.grid->v[i][j][k] /= ms.grid->m[i][j][k];
            } else {
              ms.grid->v[i][j][k] = Eigen::Matrix<float,3,1>::Zero();
            }
          }
        }
      }
      std::cout << "calculated grid velocities" << std::endl;

      // apply gravity to each node
      for (int i = 0; i < x; i++) {
        for (int j = 0; j < y; j++) {
          for (int k = 0; k < z; k++) {
            ms.grid->v[i][j][k] += gravity * dt;
          }
        }
      }
      std::cout << "applied gravity" << std::endl;

      // calculate elastic forces & velocities for grid nodes
      ms.calculateElasticForces();
      // update velocity of each node
      for (int i = 0; i < x; i++) {
        for (int j = 0; j < y; j++) {
          for (int k = 0; k < z; k++) {
            if (ms.grid->m[i][j][k] > 0.f) {
              ms.grid->v[i][j][k] += dt * ms.grid->f[i][j][k] / ms.grid->m[i][j][k];
            }
            
            // check floor and wall bounds, 5 layers
            if (i < 5 || i > x - 5) {
              ms.grid->v[i][j][k][0] = 0.f;
            }
            if (j < 5 || j > y - 5) {
              ms.grid->v[i][j][k][1] = 0.f;
            }
            if (k < 5 || k > z - 5) {
              ms.grid->v[i][j][k][2] = 0.f;
            }
          }
        }
      }
      std::cout << "applied elastic forces" << std::endl;
 
      // interpolate velocities for each particle
      for (int m = 0; m < n; m++) {
        Particle* p = ms.particles[m];
        p->v = Eigen::Matrix<float,3,1>::Zero();
        Eigen::Matrix<float,3,1> p_pos = p->x;

        int xb = int((p_pos[0] - 0.5f * h) / h);
        int yb = int((p_pos[1] - 0.5f * h) / h);
        int zb = int((p_pos[2] - 0.5f * h) / h);
        xb = std::max(0, xb);
        yb = std::max(0, yb);
        zb = std::max(0, zb);
        
        for (int i = xb; i <= xb + 2 && i < x; i++) {
          for (int j = yb; j <= yb + 2 && j < y; j++) {
            for (int k = zb; k <= zb + 2 && k < z; k++) {
              Eigen::Matrix<float,3,1> n_pos = ms.grid->getPosition(i, j, k);
              float w = ms.omega(p_pos, n_pos);
              p->v += ms.grid->v[i][j][k] * w;
            }
          }
        }
      }
      std::cout << "interpolated velocities" << std::endl;

      // update deformation gradients
      for (int m = 0; m < n; m++) {
        Particle* p = ms.particles[m];
        Eigen::Matrix<float,3,1> p_pos = p->x;

        int xb = int((p_pos[0] - 0.5f * h) / h);
        int yb = int((p_pos[1] - 0.5f * h) / h);
        int zb = int((p_pos[2] - 0.5f * h) / h);
        xb = std::max(0, xb);
        yb = std::max(0, yb);
        zb = std::max(0, zb);

        Eigen::Matrix<float,3,3> dv = Eigen::Matrix<float,3,3>::Zero();
        
        for (int i = xb; i <= xb + 2 && i < x; i++) {
          for (int j = yb; j <= yb + 2 && j < y; j++) {
            for (int k = zb; k <= zb + 2 && k < z; k++) {
              Eigen::Matrix<float,3,1> n_pos = ms.grid->getPosition(i, j, k);
              Eigen::Matrix<float,3,1> dw = ms.omegaGradient(p_pos, n_pos);
              dv += ms.grid->v[i][j][k] * dw.transpose();
            }
          }
        }

        p->F = (Eigen::Matrix<float,3,3>::Identity() + dt * dv) * p->F;
      }
      std::cout << "updated deformation gradients" << std::endl;

      // move particles
      for (int i = 0; i < n; i++) {
        Particle* p = ms.particles[i];
        p->x += p->v * dt;
      }
      std::cout << "moved particles" << std::endl;
    }
};