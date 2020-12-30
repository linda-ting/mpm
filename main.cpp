#include <Eigen/Core>
#include <Eigen/Dense>
#include <iostream>
#include <fstream>
#include <sstream>
#include <ctime>
#include <cstdlib>
#include <random>
#include <chrono>
#include <unordered_set>

#include "SimulationDriver.h"
#include "MaterialPointSystem.h"

int main(int argc, char* argv[])
{
    SimulationDriver driver;

    // set up mass spring system
    float youngs_modulus = 0.f;
    float poissons_ratio = 0.f;
    float dt = 0.f;

    // node data
    std::vector<Eigen::Matrix<float,3,1>> x;

    // segment data
    std::vector<float> rest_length;
    std::vector<Eigen::Matrix<int,3,1>> triangles;

    if (argc < 2) 
    {
        std::cout << "Please indicate test case number: 0 (cloth) or 1 (volumetric bunny)" << std::endl;
        exit(0);
    }

    if (strcmp(argv[1], "0") == 0) // two cubes case
    {
        int side = 20;
        float mult = 1.f / side;

        // set up cube A
        for (int i = 0; i <= side; i++) {
          float xc = 2.f + i * mult;
          for (int j = 0; j <= side; j++) {
            float yc = 0.5f + j * mult;
            for (int k = 0; k <= side; k++) {
              float zc = 2.f + k * mult;
              Eigen::Matrix<float,3,1> pos(xc, yc, zc);
              x.push_back(pos);
            }
          }
        }

        for (int i = 0; i < side; i++) {
          for (int j = 0; j < side; j++) {
            for (int k = 0; k < side; k++) {
              int a = (side + 1) * (side + 1) * i + (side + 1) * j + k;
              int b = (side + 1) * (side + 1) * (i + 1) + (side + 1) * j + k;
              int c = (side + 1) * (side + 1) * i + (side + 1) * j + k + 1;
              int d = (side + 1) * (side + 1) * (i + 1) + (side + 1) * j + k + 1;
              int e = (side + 1) * (side + 1) * i + (side + 1) * (j + 1) + k;
              int f = (side + 1) * (side + 1) * (i + 1) + (side + 1) * (j + 1) + k;
              int g = (side + 1) * (side + 1) * i + (side + 1) * (j + 1) + k + 1;
              int h = (side + 1) * (side + 1) * (i + 1) + (side + 1) * (j + 1) + k + 1;

              if (i == 0) {
                triangles.push_back(Eigen::Matrix<int,3,1>(a, c, e));
                triangles.push_back(Eigen::Matrix<int,3,1>(c, e, g));
              } 

              if (j == 0) {
                triangles.push_back(Eigen::Matrix<int,3,1>(a, b, c));
                triangles.push_back(Eigen::Matrix<int,3,1>(b, c, d));
              }

              if (k == 0) {
                triangles.push_back(Eigen::Matrix<int,3,1>(a, b, e));
                triangles.push_back(Eigen::Matrix<int,3,1>(b, e, f));
              }

              triangles.push_back(Eigen::Matrix<int,3,1>(b, c, e));
              triangles.push_back(Eigen::Matrix<int,3,1>(b, c, h));
              triangles.push_back(Eigen::Matrix<int,3,1>(b, d, h));
              triangles.push_back(Eigen::Matrix<int,3,1>(c, d, h));
              triangles.push_back(Eigen::Matrix<int,3,1>(c, e, h));
              triangles.push_back(Eigen::Matrix<int,3,1>(c, g, h));
              triangles.push_back(Eigen::Matrix<int,3,1>(e, g, h));
              triangles.push_back(Eigen::Matrix<int,3,1>(b, e, h));
              triangles.push_back(Eigen::Matrix<int,3,1>(b, f, h));
              triangles.push_back(Eigen::Matrix<int,3,1>(e, f, h));
            }
          }
        }

        youngs_modulus = 70.f;
        poissons_ratio = 0.4f;
        dt = 0.01;
        driver.test = "cubes";
        driver.ms.test = "cubes";
    }

    else if (strcmp(argv[1], "1") == 0) // volumetric bunny case
    { 
        std::ifstream pointsFile;
        pointsFile.open("data/points");
        std::string pointsLine;
        getline(pointsFile, pointsLine);
        int numPt = std::stoi(pointsLine.substr(0, pointsLine.find(" ")));
        for (int i = 0; i < numPt; i++) {
          getline(pointsFile, pointsLine);
          char* tok;
          char line[100];
          std::strncpy(line, pointsLine.c_str(), sizeof(line));

          tok = std::strtok(line, " ");
          std::string sx(tok);
          float xc = std::stof(sx);

          tok = std::strtok(NULL, " ");
          std::string sy(tok);
          float yc = std::stof(sy);

          tok = std::strtok(NULL, " ");
          std::string sz(tok);
          float zc = std::stof(sz);

          x.push_back(Eigen::Matrix<float,3,1>(xc, yc, zc));
        }

        std::ifstream cellsFile;
        cellsFile.open("data/cells");
        std::string cellsLine;
        getline(cellsFile, cellsLine);
        int numSeg = std::stoi(cellsLine.substr(0, cellsLine.find(" ")));
        for (int i = 0; i < numSeg; i++) {
          getline(cellsFile, cellsLine);
          char* tok;
          char line[100];
          std::strncpy(line, cellsLine.c_str(), sizeof(line));

          tok = std::strtok(line, " ");
          std::string s1(tok);
          int idx1 = std::stoi(s1);

          tok = std::strtok(NULL, " ");
          std::string s2(tok);
          int idx2 = std::stoi(s2);

          tok = std::strtok(NULL, " ");
          std::string s3(tok);
          int idx3 = std::stoi(s3);

          tok = std::strtok(NULL, " ");
          std::string s4(tok);
          int idx4 = std::stoi(s4);

          std::vector<int> idx = {idx1, idx2, idx3, idx4};
          sort(idx.begin(), idx.end());
          std::vector<Eigen::Matrix<int,3,1> >::iterator it;

          for (int j = 0; j < 4; j++) {
            for (int k = j + 1; k < 4; k++) {
              for (int l = k + 1; l < 4; l++) {
                Eigen::Matrix<int,3,1> tri(idx[j], idx[k], idx[l]);
                it = find(triangles.begin(), triangles.end(), tri);
                if (it == triangles.end()) {
                  triangles.push_back(tri);
                }
              }
            }
          }
        }

        youngs_modulus = 90.f;
        poissons_ratio = 0.35f;
        dt = 0.01;
        driver.test = "bunny";
        driver.ms.test = "bunny";
    }

    else {
        std::cout << "Wrong case number!" << std::endl;
        exit(0);
    }

    driver.dt = dt;
    driver.ms.youngs_modulus = youngs_modulus;
    driver.ms.poissons_ratio = poissons_ratio;
    driver.ms.recalculateConstants();
    driver.ms.constructMesh(x.size(), x, triangles.size(), triangles);
    driver.run(120);

    return 0;
}