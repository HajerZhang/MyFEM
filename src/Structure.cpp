#include "Structure.h"
#include <fstream>
#include <string>
#include <sstream>
#include <iostream>

using namespace std;

void tria::read_dat_file(const string& datFilePath)
{   
    ifstream datFile(datFilePath);
    if (!datFile.is_open()) {
        cerr << "Error: Unable to open dat file: " << datFilePath << endl;
        exit(1);
    }
    
    vector<vector<int>> boudary;
    struct force {
        int id;
        int direct;
        double value;
    };
    vector<force> fo;

    string line, subline;
    while (getline(datFile, line)) {
        if (line.empty()) {
            continue;
        }

        istringstream iss(line);
        string key;
        iss >> key;

        if (key == "nodes:") {
            iss >> node_sum;
            ndof = node_sum * 2;
            point.resize(node_sum);
            for (int i = 0; i < node_sum; i++) {
                nodes point[i];
            }
        } else if (key == "nodes_list:") {
            for (int i = 0; i < node_sum; ++i) {
                getline(datFile, subline);
                istringstream is(subline);
                is >> point[i].id >> point[i].coordinate[0] >> point[i].coordinate[1];
            }
        } else if (key == "elements:") {
            iss >> cell_sum;
            cell.resize(cell_sum);
            for (int i = 0; i < cell_sum; i++) {
                cell_structure cell[i];
            }
        } else if(key == "elements_list:") {
            for (int i = 0; i < cell_sum; ++i) {
                getline(datFile, subline);
                istringstream is(subline);
                is >> cell[i].id >> cell[i].cell_node[0] >> cell[i].cell_node[1] >> cell[i].cell_node[2] >> cell[i].cell_node[3];
            }
        } else if (key == "E:") {
            iss >> E;
        } else if (key == "nu:") {
            iss >> nu;
        } else if (key == "gauss:") {
            iss >> gauss_num;
        } else if (key == "boundary:") {
            while (getline(datFile, line)) {
                istringstream is(line);
                int value;
                if (!(is >> value)) {
                    datFile.seekg(-(static_cast<streamoff>(line.length()) + 1), ios::cur);
                break;
                } else {
                    istringstream js(line);
                    boudary.push_back(vector<int>(2));
                    js >> boudary.back()[0] >> boudary.back()[1];
                }
            }
        } else if (key == "force:" ){
            while (getline(datFile, line)) {
                istringstream is(line);
                int value;
                if (!(is >> value)) {
                    datFile.seekg(-(static_cast<streamoff>(line.length()) + 1), ios::cur);
                break;
                } else {
                    istringstream js(line);
                    fo.push_back(force());
                    js >> fo.back().id >> fo.back().direct >> fo.back().value;
                }
            }
        } else if (key == "end") {
            break;
        } 
    }
    datFile.close();

    // // output for data check
    // cout << "node_sum: " << node_sum << endl;
    // cout << "cell_sum: " << cell_sum << endl;
    // cout << "E: " << E << endl;
    // cout << "nu: " << nu << endl;
    // cout << "gauss_num: " << gauss_num << endl;
    // for (int i = 0; i < node_sum; ++i) {
    //     cout << point[i].id << " " << point[i].coordinate[0] << " " << point[i].coordinate[1] << endl;
    // }
    // for (int i = 0; i < cell_sum; ++i) {
    //     cout << cell[i].id << " " << cell[i].cell_node[0] << " " << cell[i].cell_node[1] << " " << cell[i].cell_node[2] << " " << cell[i].cell_node[3] << endl;
    // }
    // for (int i = 0; i < boudary.size(); ++i) {
    //     cout << boudary[i][0] << " " << boudary[i][1] << endl;
    // }
    // for (int i = 0; i < fo.size(); ++i) {
    //     cout << fo[i].id << " " << fo[i].direct << " " << fo[i].value << endl;
    // }
    
    // element initial
    for (int i = 0; i < cell_sum; i++) {
        cell[i].Initial(i, E, nu, gauss_num, point, cell[i].cell_node);
    }

    // load setting
    F.resize(ndof);
    for (int i = 0; i < fo.size(); ++i) {
        F.coeffRef(-1 + fo[i].id * 2 - (2 - fo[i].direct)) = fo[i].value;
    }

    // boundary condition
    Eigen::VectorXd xx = Eigen::VectorXd::Ones(ndof);
	Eigen::VectorXd yy = Eigen::VectorXd::Zero(ndof);
    for (int i = 0; i < boudary.size(); ++i) {
        xx.coeffRef(-1 + boudary[i][0] * 2 - (2 - boudary[i][1])) = 0;
        yy.coeffRef(-1 + boudary[i][0] * 2 - (2 - boudary[i][1])) = 1;
    }
    diagonalmatrix_one = xx.asDiagonal();
	diagonalmatrix_two = yy.asDiagonal();

}

void cell_structure::Initial(int& it, double& E, double& nu, int& gaus_sum,
    const vector<nodes>& node, const vector<int>& cell)
{
    
    
}

nodes::nodes()
{
    coordinate.resize(2);
}

nodes::~nodes()
{

}

cell_structure::cell_structure()
{
    cell_node.resize(4);
}

cell_structure::~cell_structure()
{

}

tria::tria()
{

}

tria::~tria()
{

}