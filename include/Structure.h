#pragma once
#include <vector>
#include <memory>
using namespace std;

class nodes
{
public:
    int nodes_id;
    vector<double> node_coordinate;
       
};

class cell_face // means line in 2D
{
public:
    int face_num; // means line number of the face
    vector<int> face;
};

class cell_structure
{
public:
    vector<unique_ptr<cell_face>> face;
};

