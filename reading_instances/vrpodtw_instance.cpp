/******************************************************************************
 * vrpodtw_instance.cpp: implementation for VRPODTW_Instance class.
 *****************************************************************************/

#include "vrpodtw_instance.hpp"
#include <fstream>
#include <stdexcept>
#include <iostream>

using namespace std;

//-----------------------------[ Constructor ]--------------------------------//

VRPODTW_Instance::VRPODTW_Instance(const std::string& filename):
    num_nodes(0),
    num_drivers(0),
    num_od(0),
    capacities(),
    TW_matrix(),
    initial_dist(),
    demands(),
    distances()
{
    ifstream file(filename, ios::in);
    if(!file)
        throw runtime_error("Cannot open instance file");

    file.exceptions(ifstream::failbit | ifstream::badbit);
    try {
        file >> num_nodes;
        file >> num_drivers;
        file >> num_od;

        capacities.resize(num_drivers);
        for(unsigned i = 0; i < num_drivers; ++i)
            file >> capacities[i];

        distances.resize(((num_nodes-num_drivers+num_od) * (num_nodes-num_drivers+num_od + 1)) / 2-num_od * (num_od+1)/2);
        for(unsigned i = 0; i < distances.size(); ++i)
            file >> distances[i];

        TW_matrix.resize(num_nodes);
        for(unsigned i = 0; i < num_nodes; ++i)
            TW_matrix[i].resize(2);
        for(int i = 0; i < num_nodes; i++){
            file >> TW_matrix[i][0];
            file >> TW_matrix[i][1];}

        initial_dist.resize(num_drivers);
        for(unsigned i = 0; i < num_drivers; ++i)
            file >> initial_dist[i];

        demands.resize(num_nodes-num_drivers);
        for(unsigned i = 0; i < demands.size(); ++i)
            file >> demands[i];
    }
    catch(std::ifstream::failure& e) {
        throw fstream::failure("Error reading the instance file");
    }
}
//-------------------------------[ Distance ]---------------------------------//

double VRPODTW_Instance::distance(int i, int j) const {
    if(i > j)
        swap(i, j);
    return distances[(i * (num_nodes-num_drivers+num_od)) - ((i-1) * i / 2) +
            (j-i-1)];
}
