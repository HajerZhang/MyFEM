#include "Structure.h"
#include "gauss.h"

using namespace std;

void tria::Read_dat_file(const string& datFilePath)
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
    //     cout << point[i].id << " " << point[i].coordinate[0] << " " << point[i].coordinate[1] << " 0"<< endl;
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
        cell[i].Initial(E, nu, gauss_num, point, cell[i].cell_node);
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


void cell_structure::Initial(double& E, double& nu, int& gaus_sum,
    const vector<nodes>& node, const vector<int>& cell)
{
    ke.resize(pow(gaus_sum, 2), vector<vector<double>>(8, vector<double>(8, 0.0)));
    reflect.resize(8, 0);
    location.resize(pow(gaus_sum, 2), vector<double>(2, 0.0));
    gauss_iso* aaa;
    vector<double> x1{ node[cell[0] - 1].coordinate[0], node[cell[1] - 1].coordinate[0], 
                        node[cell[2] - 1].coordinate[0], node[cell[3] - 1].coordinate[0] };
	vector<double> y1{ node[cell[0] - 1].coordinate[1], node[cell[1] - 1].coordinate[1], 
                        node[cell[2] - 1].coordinate[1], node[cell[3] - 1].coordinate[1] };
	vector<double> z1;
    aaa = new gauss_iso(ke, E, nu, gaus_sum, x1, y1, z1);
    aaa->compute(location);
    vector<vector<double>> keke(8, vector<double>(8, 0.0));
	for (int i = 0; i < 8; i++) {
		for (int j = 0; j < 8; j++) {
			keke[i][j] += ke[0][i][j] + ke[1][i][j] + ke[2][i][j] + ke[3][i][j];
		}
	}
    reflect[0] = 2 * cell[1] - 1;
	reflect[1] = 2 * cell[1];
	reflect[2] = 2 * cell[2] - 1;
	reflect[3] = 2 * cell[2];
	reflect[4] = 2 * cell[3] - 1;
	reflect[5] = 2 * cell[3];
	reflect[6] = 2 * cell[0] - 1;
	reflect[7] = 2 * cell[0];
}

void tria::ComputeObjective()
{
    int gaussnum = pow(2, 2);  //高斯积分点数
    int cellndofs = 2 * pow(2, 2);  //一个单元的自由度数
    typedef Eigen::Triplet<double> T;
    vector<T> tripletList;
    tripletList.reserve(cell_sum * gaussnum * cellndofs * cellndofs);
    for (int itt = 0; itt < cell_sum; itt++) {
        //遍历itt号单元的单刚所有元素
        for (int k = 0; k < gaussnum; k++) {
            for (int i = 0; i < cellndofs; i++) {
                for (int j = 0; j < cellndofs; j++) {
                    tripletList.push_back(T(cell[itt].reflect[i] - 1, cell[itt].reflect[j] - 1,   E * cell[itt].ke[k][i][j]));
                }
            }
        }
    }
    Eigen::SparseMatrix<double> K_global(ndof, ndof);

    K_global.reserve(ndof * 81);
    K_global.setFromTriplets(tripletList.begin(), tripletList.end());
    K_global = diagonalmatrix_one * K_global * diagonalmatrix_one + diagonalmatrix_two;

    Eigen::SparseVector<double> U;

    Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>> solver;
    solver.analyzePattern(K_global);
    solver.factorize(K_global);
    U = solver.solve(F);

    for(int i = 0; i < node_sum; i++)
    {
        point[i].displacement[0] = U.coeffRef(2 * i);
        point[i].displacement[1] = U.coeffRef(2 * i + 1);
    }

}

void tria::Output()
{
    using namespace std;
    
    const char* folderPath = "./Results";

    // Check if the "data" folder exists, create it if not
    struct stat info;
    if (stat(folderPath, &info) != 0 || !(info.st_mode & S_IFDIR)) {
        // Folder doesn't exist, create it
        if (mkdir(folderPath, 0777) != 0) {
            cerr << "Error creating directory: " << folderPath << std::endl;
            return;
        }
    }

    string filename = "./Results/result.vtk";
    ofstream outputFile(filename);

    if (!outputFile.is_open()) {
        cerr << "Error opening file: " << filename << std::endl;
        return;
    }

    outputFile << "# vtk DataFile Version 3.0" << endl;
    outputFile << "vtk output" << endl;
    outputFile << "ASCII" << endl;
    outputFile << "DATASET POLYDATA" << endl;
    outputFile << "POINTS " << node_sum << " float" << endl;
    for(int i = 0; i < node_sum; i++)
    {
        outputFile << point[i].coordinate[0] << " " << point[i].coordinate[1] << " " << 0 << endl;
    }
    outputFile << "POLYGONS " << cell_sum << " " << cell_sum * 5 << endl;
    for(int i = 0; i < cell_sum; i++)
    {   
        // number start from 0 in vtk file, need to minus 1
        outputFile << 4 << " " << cell[i].cell_node[0]-1 << " " << cell[i].cell_node[1]-1 << " " << cell[i].cell_node[2]-1 << " " << cell[i].cell_node[3]-1 << endl;
    }
    outputFile << "POINT_DATA " << node_sum << endl;
    outputFile << "VECTORS displacement float" << endl;
    for(int i = 0; i < node_sum; i++)
    {
        outputFile << point[i].displacement[0] << " " << point[i].displacement[1] << " " << 0 << endl;
    }
    // outputFile << "CELL_DATA " << cell_sum << endl;
    // outputFile << "SCALARS density float" << endl;
    // outputFile << "LOOKUP_TABLE default" << endl;
    // for(int i = 0; i < nele; i++)
    // {
    //     outputFile << x(i) << endl;
    // }
    outputFile.close();
}

nodes::nodes()
{
    coordinate.resize(2);
    displacement.resize(2);
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