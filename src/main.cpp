#include "Structure.h"
#include <iostream>

void deletePrevious(const std::string& directoryPath);

int main()
{   
    using namespace std;
    // data folder path
    string directoryPath = "Results";
    string inputdataPath = "./Input/";
    // delete previous vtk files
    deletePrevious(directoryPath);
    // structure object
    tria triangulation;
    // read dat file
    triangulation.Read_dat_file(inputdataPath+"test.dat");
    // compute of linear elastic structure
    triangulation.ComputeObjective();
    triangulation.Output();
    return 0;
}

void deletePrevious(const std::string& directoryPath) {
    DIR* dir = opendir(directoryPath.c_str());
    if (dir) {
        struct dirent* entry;
        while ((entry = readdir(dir)) != nullptr) {
            std::string fileName = entry->d_name;
            if (fileName.length() > 4 &&
                fileName.substr(0, 7) == "result_" &&
                fileName.substr(fileName.length() - 4) == ".vtk") {
                std::string filePath = directoryPath + "/" + fileName;
                std::remove(filePath.c_str());
            }
        }
        closedir(dir);
    }
}