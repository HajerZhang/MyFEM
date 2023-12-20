#include "Structure.h"
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>


int main()
{   
    using namespace std;
    int num_of_nodes;
    string datFilePath = "mesh_model.dat";
    vector<unique_ptr<nodes>> point;
    ifstream datFile(datFilePath);
    if (!datFile.is_open()) {
        cerr << "Error: Unable to open dat file: " << datFilePath << endl;
        return 1;
    }

    string line, subline;
    while (getline(datFile, line)) {
        if (line.empty()) {
            continue;
        }

        istringstream iss(line);
        string key;
        iss >> key;

        if (key == "nodes:") {
            iss >> num_of_nodes;
            point.reserve(num_of_nodes);
        } else if (key == "nodes_list:") {
             for (int i = 0; i < num_of_nodes; ++i) {
                getline(datFile, subline);
                istringstream is(subline);
                point.push_back(make_unique<nodes>());
                is >> point[i]->nodes_id >> point[i]->node_coordinate[0] >> point[i]->node_coordinate[1];
            }
        } 
    }
    for (int i = 0; i < num_of_nodes; ++i) {
        cout << point[i]->nodes_id << " " << point[i]->node_coordinate[0] << " " << point[i]->node_coordinate[1] << endl;
    }
    return 0;
}
