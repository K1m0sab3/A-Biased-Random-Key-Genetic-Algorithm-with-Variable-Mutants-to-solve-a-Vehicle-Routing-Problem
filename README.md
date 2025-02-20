# A-Biased-Random-Key-Genetic-Algorithm-with-Variable-Mutants-to-solve-a-Vehicle-Routing-Problem
BRKGA decoder, configuration files, instance reader, dataset and main file for the VRPODTW.

This repository contains the decoder class for the BRKGA, parameter configuration files, instance reading procedure, dataset and the main execution file for
solving the Vehicle Routing Problem with Occasional Drivers (VRPOD) studied in Festa et al. “A Biased Random-Key Genetic Algorithm with Variable Mutants to solve
a Vehicle Routing Problem”, under review.

- In “**Instances**” folder the four sets of instances used in the paper are provided, along with a README file presenting the format and characteristics  of the
  instances.
  
- In the “**reading_instances**” folder, the implementation for reading problem instances is provided. The class VRPODTW_Instance is responsible for parsing instance
  files, extracting key problem parameters, and storing relevant data structures. It reads information such as the number of nodes, drivers, and occasional
  drivers, vehicle capacities, time windows, initial distances, demands, and the distance matrix. The class ensures efficient data handling and provides a
  function to retrieve the distance between any two nodes.
  
- In “**decoder**” folder there is the VRPODTW_Decoder class, which implements the BRKGA decoding process for the problem. In particular, it assigns customers to
  Occasional Drivers (OD) or company couriers while ensuring feasibility based on time windows and capacity constraints. The decoder processes customers in a
  sorted order based on chromosome values, attempting to insert them into an OD path while minimizing cost. If no feasible insertion is possible, the customer is
  assigned to company couriers.
  
- In “**configurations**” folder there are the parameter settings used for running the BRKGA algorithm. There are three configuration files.
  >  config_Test&Small: parameters tuned for small-sized instances;
  > config_Medium: parameters tuned for medium-sized instances;
  > config_Large: parameters tuned for large-sized instances.
   Each file defines key evolutionary parameters such as population size, elite and mutant percentages, parent selection strategy, and
  path-relinking settings.
  These configurations ensure the BRKGA is properly calibrated for different problem scales.
  
- “**main.cpp**” is the file that integrates all components and runs the BRKGA for solving VRPOD instances.
  > It loads an instance file and parses it using VRPODTW_Instance.
  > It reads BRKGA configuration parameters and initializes the algorithm.
  > It evolves the population while applying mutation, selection, and path relinking strategies.
  > It monitors stopping criteria (generations and improvement).
  > It outputs key performance metrics, including cost, iteration count, and runtime statistics.
  The main file is designed to run multiple instances with different seeds and configurations, making it suitable for computational
  experiments and performance evaluation.
  
**Implementation Details ** 
The classes are implemented in C++ and are compatible with the BRKGA-MP-IPR library, available at:  
[BRKGA-MP-IPR (C++ version)](https://github.com/ceandrade/brkga_mp_ipr_cpp)


**Usage**  
Integrate the provided decoder into the BRKGA framework and upload the configuration files accordingly. The instance reader allows processing problem instances in the adapted format.

For any inquiries, you can contact the author Edoardo Scalzo at edoardo.scalzo@unical.it.
