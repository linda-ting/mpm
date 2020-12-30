#pragma once

#include <Eigen/Core>
#include <Eigen/Dense>
#include <vector>
#include <string>
#include <cstdlib>

class Particle {
public:
    float m;
    Eigen::Matrix<float,3,1> x;
    Eigen::Matrix<float,3,1> v;
    Eigen::Matrix<float,3,3> F;
    Eigen::Matrix<float,3,3> P;

    Particle(float m, Eigen::Matrix<float,3,1> x) {
        this->m = m;
        this->x = x;
        this->v = Eigen::Matrix<float,3,1>::Zero();
        this->F = Eigen::Matrix<float,3,3>::Identity();
        this->P = Eigen::Matrix<float,3,3>::Zero();
    }
};