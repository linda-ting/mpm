#pragma once

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/SVD>

#include <vector>
#include <string>
#include <fstream>
#include <cstdlib>

#include "grid.h"
#include "particle.h"
#include "mesh_query.h"

class MaterialPointSystem{
public:
    Grid* grid;
    std::vector<Particle*> particles;
    int n;
    MeshObject* mesh;
    float youngs_modulus;
    float poissons_ratio;
    float mu;
    float lambda;
    float particle_volume;
    std::string test;

    MaterialPointSystem(int n) {
      this->grid = new Grid();
      this->n = n;
      this->particle_volume = 1.f / n;
    }

    void recalculateConstants() {
      this->mu = youngs_modulus / (1.f + poissons_ratio);
      this->lambda = youngs_modulus * poissons_ratio / ((1.f + poissons_ratio) * (1.f - 2.f * poissons_ratio));
    }

    void constructMesh(int num_vert, std::vector<Eigen::Matrix<float,3,1>> vertices,
                       int num_tri, std::vector<Eigen::Matrix<int,3,1>> triangles) {
      double pos[3 * num_vert];
      for (int i = 0; i < num_vert; i++) {
        Eigen::Matrix<float,3,1> p = vertices[i];
        pos[3 * i] = double(p[0]);
        pos[3 * i + 1] = double(p[1]);
        pos[3 * i + 2] = double(p[2]);
      }

      int tri[3 * num_tri];
      for (int i = 0; i < num_tri; i++) {
        Eigen::Matrix<int,3,1> t = triangles[i];
        tri[3 * i] = double(t[0]);
        tri[3 * i + 1] = double(t[1]);
        tri[3 * i + 2] = double(t[2]);
      }

      this->mesh = construct_mesh_object(num_vert, pos, num_tri, tri);
      sampleParticles();
      std::cout << "done sampling particles" << std::endl;
    }

    void sampleParticles() {
      int numParticle = 0;
      Eigen::Matrix<float,3,1> o = grid->origin;
      float xlim = grid->dims[0] * grid->h;
      float ylim = grid->dims[1] * grid->h;
      float zlim = grid->dims[2] * grid->h;
      std::ofstream obj;
      std::string filename = "data/" + test + "_" + std::to_string(n) + ".obj";
      obj.open(filename, std::ios::trunc | std::ios::binary);

      while (numParticle < n) {
        float x = xlim * std::rand() / RAND_MAX;
        float y = ylim * std::rand() / RAND_MAX;
        float z = zlim * std::rand() / RAND_MAX;
        float mass = 10.f / n;
        Eigen::Matrix<float,3,1> pos(x, y, z);

        double d_pos[3];
        d_pos[0] = pos[0];
        d_pos[1] = pos[1];
        d_pos[2] = pos[2];
        if (point_inside_mesh(d_pos, mesh)) {
          Particle* p = new Particle(mass, pos);
          particles.push_back(p);
          numParticle++;

          std::string v = "v " + std::to_string(x) + " " + std::to_string(y) + " " + std::to_string(z) + "\n";
          obj << v;
        }
      }
      obj.close();
    }

    float nHat(float p, float n) {
      float x = abs((p - n) / grid->h);
      if (x >= 0 && x < 0.5f) {
        return 0.75f - pow(x, 2);
      } else if (x >= 0.5f && x < 1.5f) {
        return 0.5f * pow(1.5f - x, 2);
      }
      return 0.f;
    }

    float nHatPrime(float p, float n) {
      float x = (p - n) / grid->h;
      float abs_x = abs(x);
      if (abs_x >= 0 && abs_x < 0.5f) {
        return -2.f * x;
      } else if (abs_x >= 0.5f && abs_x < 1.5f) {
        return x * (1.f - 3.f / (2.f * abs_x));
      }
      return 0.f;
    }

    float omega(Eigen::Matrix<float,3,1> p_pos, Eigen::Matrix<float,3,1> n_pos) {
      float nx = nHat(p_pos[0], n_pos[0]);
      float ny = nHat(p_pos[1], n_pos[1]);
      float nz = nHat(p_pos[2], n_pos[2]);
      return nx * ny * nz;
    }

    Eigen::Matrix<float,3,1> omegaGradient(Eigen::Matrix<float,3,1> p_pos, Eigen::Matrix<float,3,1> n_pos) {
      Eigen::Matrix<float,3,1> dw;
      dw[0] = nHatPrime(p_pos[0], n_pos[0]) * nHat(p_pos[1], n_pos[1]) * nHat(p_pos[2], n_pos[2]) / grid->h;
      dw[1] = nHat(p_pos[0], n_pos[0]) * nHatPrime(p_pos[1], n_pos[1]) * nHat(p_pos[2], n_pos[2]) / grid->h;
      dw[2] = nHat(p_pos[0], n_pos[0]) * nHat(p_pos[1], n_pos[1]) * nHatPrime(p_pos[2], n_pos[2]) / grid->h;
      return dw;
    }

    void calculateStress(Particle* p) {
      Eigen::JacobiSVD<Eigen::Matrix<float,3,3>> svd(p->F, Eigen::ComputeFullV | Eigen::ComputeFullU);
      Eigen::Matrix<float,3,3> U = svd.matrixU();
      Eigen::Matrix<float,3,3> V = svd.matrixV();
      Eigen::Matrix<float,3,3> R = U * V.transpose();
      float J = p->F.determinant();
      Eigen::Matrix<float,3,3> P = 2 * mu * (p->F - R) +
                                   lambda * (J - 1) * J * p->F.inverse().transpose();
      p->P = P;
    }

    void calculateElasticForces() {
      int x = grid->dims[0];
      int y = grid->dims[1];
      int z = grid->dims[2];
      float h = grid->h;

      // for all particles
      for (int m = 0; m < n; m++) {
        Particle* p = particles[m];
        Eigen::Matrix<float,3,1> p_pos = p->x;              
        calculateStress(p);

        int xb = int((p_pos[0] - 0.5f * h) / h);
        int yb = int((p_pos[1] - 0.5f * h) / h);
        int zb = int((p_pos[2] - 0.5f * h) / h);
        xb = std::max(0, xb);
        yb = std::max(0, yb);
        zb = std::max(0, zb);
        
        // for all nearby grid nodes, add corresponding force
        for (int i = xb; i <= xb + 2 && i < x; i++) {
          for (int j = yb; j <= yb + 2 && j < y; j++) {
            for (int k = zb; k <= zb + 2 && k < z; k++) {
              Eigen::Matrix<float,3,1> n_pos = grid->getPosition(i, j, k);
              Eigen::Matrix<float,3,1> dw = omegaGradient(p_pos, n_pos);      
              grid->f[i][j][k] -= particle_volume * p->P * p->F.transpose() * dw;
            }
          }
        }
      }
    }

    void dumpPoly(std::string filename)
    {
        std::ofstream fs;
        fs.open(filename);
        fs << "POINTS\n";
        int count = 0;
        for (auto P : particles) {
            Eigen::Matrix<float,3,1> X = P->x;
            fs << ++count << ":";
            for (int i = 0; i < 3; i++)
                fs << " " << (double)X(i);
            fs << "\n";
        }
        fs << "POLYS\n";
        fs << "END\n";
        fs.close();
    }
};