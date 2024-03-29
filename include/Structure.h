#pragma once
#include <vector>
#include <memory>
#include <cmath>

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <Eigen/IterativeLinearSolvers>

#include <string>
#include <dirent.h>
#include <sys/stat.h>
#include <iostream>
#include <fstream>

using namespace std;

class nodes
{
public:
    int id;
    vector<double> coordinate;
    vector<double> displacement;
    nodes();
    ~nodes();
};

// class cell_face // means line in 2D
// {
// public:
//     int face_num; // means line number of the face
//     vector<int> face_node;

// };

class cell_structure
{
public: 
    int id;
    vector<int> cell_node;
    vector<vector<vector<double>>> ke;
    // vector<cell_face> face;
    vector<int> reflect;

    cell_structure();
    ~cell_structure();
    void Initial(double& E, double& nu, int& gaus_sum,
    const vector<nodes>& node, const vector<int>& cell);

    vector<vector<double>> location;
};

class tria
{
public:
    tria();
    ~tria();

    double E;
    double nu;
    int gauss_num;

    int cell_sum, node_sum;

    int ndof;

    vector<nodes> point;
    vector<cell_structure> cell;

    Eigen::SparseVector<double> F;
    
    Eigen::SparseMatrix<double> diagonalmatrix_one;
    Eigen::SparseMatrix<double> diagonalmatrix_two;

    void Read_dat_file(const string& datFilePath);
    void Compute();
    void Output();
};

