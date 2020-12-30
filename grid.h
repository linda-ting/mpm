#pragma once

#include <Eigen/Core>
#include <Eigen/Dense>
#include <vector>
#include <string>
#include <cstdlib>

class Grid {
public:
    Eigen::Matrix<float,3,1> origin;
    float h;
    Eigen::Matrix<float,3,1> dims;
    std::vector<std::vector<std::vector<float>>> m;
    std::vector<std::vector<std::vector<Eigen::Matrix<float,3,1>>>> v;
    std::vector<std::vector<std::vector<Eigen::Matrix<float,3,1>>>> f;

    Grid() {
      // grid goes from (0,0,0) to (5,5,5)
      this->origin = Eigen::Matrix<float,3,1>(0, 0, 0);
      this->h = 0.05f;
      this->dims = Eigen::Matrix<float,3,1>(101, 101, 101);

      int x = dims[0];
      int y = dims[1];
      int z = dims[2];
      m = std::vector<std::vector<std::vector<float>>>(x);
      v = std::vector<std::vector<std::vector<Eigen::Matrix<float,3,1>>>>(x);
      f = std::vector<std::vector<std::vector<Eigen::Matrix<float,3,1>>>>(x);
      for (int i = 0; i < x; i++) {
        m[i] = std::vector<std::vector<float>>(y);
        v[i] = std::vector<std::vector<Eigen::Matrix<float,3,1>>>(y);
        f[i] = std::vector<std::vector<Eigen::Matrix<float,3,1>>>(y);
        for (int j = 0; j < y; j++) {
          m[i][j] = std::vector<float>(z);
          v[i][j] = std::vector<Eigen::Matrix<float,3,1>>(z);
          f[i][j] = std::vector<Eigen::Matrix<float,3,1>>(z);
          for (int k = 0; k < z; k++) {
            m[i][j][k] = 0.f;
            v[i][j][k] = Eigen::Matrix<float,3,1>::Zero();
            f[i][j][k] = Eigen::Matrix<float,3,1>::Zero();
          }
        }
      }
    }

    Eigen::Matrix<float,3,1> getPosition(int xi, int yi, int zi) {
      return Eigen::Matrix<float,3,1>(xi * h, yi * h, zi * h);
    }
};