/******************************************************************************
 * vprodtw_instance.hpp: interface for VRPODTW_Instance class.
 *****************************************************************************/

#ifndef VRPODTW_INSTANCE_HPP_
#define VRPODTW_INSTANCE_HPP_

#include <string>
#include <vector>

class VRPODTW_Instance {
public:
    /// Default Constructor.
    VRPODTW_Instance(const std::string& filename);

    /// Return the distance between nodes i and j.
    double distance(int i, int j) const;

public:
    /// Number of elements.
    unsigned num_nodes;
    unsigned num_drivers;
    unsigned num_od;

    /// Vector of capacities.
    std::vector<int> capacities;

    /// Distances between the nodes.
    std::vector<double> distances;

    /// Matrix of time windows of customers and drivers
    std::vector<std::vector<int>> TW_matrix;

    /// Distances between od origin and central depot
    std::vector<double> initial_dist;

    /// Vector of customer demands
    std::vector<int> demands;

};

#endif // VRPODTW_INSTANCE_HPP_
