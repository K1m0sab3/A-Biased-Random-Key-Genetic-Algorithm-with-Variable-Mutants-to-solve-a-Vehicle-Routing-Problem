/***********************************************************************************
 * vrpodtw_decoder.cpp: inserts the processed customer at the end of the Occasional
 * Driver (OD) path, provided that constraints such as time windows and capacity allow it.
 * If the customer cannot be inserted into the current OD path due to these constraints,
 * the algorithm attempts to insert them into the next OD path. If the customer cannot
 * be inserted in any OD path, they are added to the list of customers to be served
 * by company couriers.
 *****************************************************************************/

#include "vrpodtw_decoder.hpp"
#include <algorithm>
#include <vector>
#include "iostream"

using namespace std;
using namespace BRKGA;

//-----------------------------[ Constructor ]--------------------------------//

VRPODTW_Decoder::VRPODTW_Decoder(const VRPODTW_Instance& _instance, int _pr_delivery_value):
    instance(_instance),
    pr_delivery_value(_pr_delivery_value),
    paths()
{}

//-------------------------------[ Decode ]-----------------------------------//

double VRPODTW_Decoder::decode(Chromosome& chromosome, bool /* not-used */, double currentBestFitness) {
    //initializing data structures
    double cost = 0;
    int num_drivers = instance.num_drivers;
    int num_customers = instance.num_nodes-instance.num_drivers;
    int customers_already_processed = num_customers;
    int num_occasional_drivers = instance.num_od;
    vector<int> capacity = instance.capacities;
    vector<vector<int>> TWs = instance.TW_matrix;
    vector<double> initial_dist = instance.initial_dist;
    const vector<int>& demands= instance.demands;
    vector<int> initial_capacity = instance.capacities;
    vector<vector<int>> local_paths(num_drivers);
    paths = local_paths;

    for(int i = 0; i < num_drivers; i++)//initialization driver paths
        paths[i].push_back(num_customers);

    vector<pair<double, unsigned>> customer_processed(num_customers);//for permutation of customers
    for(int i = 0; i < num_customers; i++)
        customer_processed[i] = make_pair(chromosome[i], i);
    sort(customer_processed.begin(), customer_processed.end());//make a customer permutation

    vector<pair<double, unsigned>> driver_processed(num_drivers);//for permutation of drivers
    for(int i = 0; i < num_drivers; i++)
        driver_processed[i] = make_pair(chromosome[i+num_customers], i);
    sort(driver_processed.begin(), driver_processed.end());//make a driver permutation

    double pr_delivery, d_delivery;
    for(int i = 0; i < num_customers; i++){//start processing phase for customers
        int customer = customer_processed[i].second;

        for(int j = 0; j < num_drivers; j++){
            int driver = driver_processed[j].second;
            int dist = int(instance.distance(customer, paths[driver][paths[driver].size()-1]));
            pr_delivery = rand() % pr_delivery_value;//assigns a probability threshold for delivery selection
            if (capacity[driver] >= (demands[customer]) && pr_delivery != 0 &&
                    TWs[num_customers+driver][1] >= TWs[customer][0] &&
                            TWs[num_customers+driver][0] + dist <= TWs[customer][1] &&
                            TWs[num_customers+driver][0] + dist <= TWs[num_customers+driver][1]) {//constraint control
                //data update
                capacity[driver] -= demands[customer];
                customers_already_processed--;
                if (paths[driver].size() == 1) cost -= int(initial_dist[driver]);
                TWs[num_customers + driver][0] = max(TWs[num_customers + driver][0] + dist, TWs[customer][0]);
                paths[driver].push_back(customer);
                if (driver >= num_occasional_drivers)
                    cost += dist;
                else
                    cost += dist * COMPENSATION;
                break;//customer processed and served, break the cycles
            }
        }
    }
    if (customers_already_processed > 0)//force the chromosome to be ejected from evolution if there is at least one unserved customer
        return INT_MAX;

    //add the final distance, from the last customer served to the destination node
    for (int od = 0; od <  num_occasional_drivers; od++){
        if (paths[od].size() != 1)
            cost += int(instance.distance(num_customers+od+1, paths[od][paths[od].size()-1])) * COMPENSATION;
    }
    for (int cd = num_occasional_drivers; cd < num_drivers; cd++){
        if (paths[cd].size() != 1)
            cost +=  int(instance.distance(paths[cd][paths[cd].size()-1], num_customers));
    }
