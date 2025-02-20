/******************************************************************************
 * main.cpp: code for calling BRKGA algorithms to solve instances of the VRPODTW.
 *****************************************************************************/

#include "reading_instances/vrpodtw_instance.hpp"
#include "decoders/vrpodtw_decoder.hpp"
#include "../../brkga_mp_ipr/brkga_mp_ipr.hpp"

#include <typeinfo>
#include <chrono>
#include <ctime>
#include <iomanip>
#include <memory>
#include <stdexcept>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
//*********************************************

using namespace std;
using namespace chrono;

//-------------------------[ Some control constants ]-------------------------//

//controls stop criteria.
enum class StopRule {
    GENERATIONS = 'G',
    IMPROVEMENT = 'I',
    TARGET = 'T',
    UNKNOWN = 'U',
};

//-------------------------[ Simple logging function ]------------------------//

void log(const string& message) {
    auto start_time = system_clock::to_time_t(system_clock::now());
    string ss(ctime(&start_time));
    ss.pop_back();  //workaround to skip unwanted end-of-line.
    cout << "\n[" << ss << "] " << message << endl;
}

//-------------------------[ Compute the elapse time ]------------------------//

double elapsedFrom(const time_point<system_clock>& start_time) {
    return chrono::duration_cast<milliseconds>(
        system_clock::now() - start_time
    ).count() / 1000.0;
}

//-------------------------------[ Main ]------------------------------------//

int main(int argc, char *argv[]) {
    const std::string USAGE =
            R"(Usage:
  main -a <instance_file> -c <seed> -d <num_run>
 )";///* <instance_file>: The path to the instance file (without the .dat extension);
    ///* <seed>: The random seed for reproducibility;
    ///* <num_run>: The identifier for the current execution run (useful for batch processing and logging).

    //declaring and initializing main variables
    string config_file;
    unsigned seed = 0;
    StopRule stop_rule = StopRule::UNKNOWN;
    double stop_arg = 0.0;
    double max_time = 0.0;
    string instance_file;
    unsigned num_threads = 8;
    bool perform_evolution = true;
    char* p;
    int num_run;

    config_file = "../configurations/config_Test&Small.conf";
    instance_file = string(argv[1]) + ".dat";;
    seed =  strtol(argv[2], &p, 10);
    num_run = atof(argv[3]);

    //setting stop rules and initial parameters
    stop_rule = static_cast<StopRule>('G');
    stop_arg = 1000;
    max_time = 900;

    int population_size = 0;
    double elite_percentage = 0;
    double mutants_percentage = 0;
    double mutants_increase_percentage = 0;
    int num_elite_parents = 0;
    int total_parents = 0;
    int num_independent_populations;
    int pr_number_pairs;
    double pr_minimum_distance;
    double pr_percentage;
    int pr_delivery;
    int num_exchange_indivuduals;
    int M_interval_between_2_improvs = 0;

    try {

        if(stop_rule != StopRule::GENERATIONS &&
           stop_rule != StopRule::IMPROVEMENT &&
           stop_rule != StopRule::TARGET)
            throw std::logic_error("Incorrect stop rule. Must be either "
                                   "(G)enerations, (I)terations, or (T)arget");

        stringstream ss;
        if(stop_rule != StopRule::TARGET && stop_arg < 1e-6) {
            ss << "when using either 'Generations' or 'Improvement' rules, "
               << "stop_arg must be positive number. Given '"
               << stop_arg << "'";
            throw std::logic_error(ss.str());
        }

        if(max_time < 1e-6) {
            ss << "max_time must be > 0. Given '" << max_time << "'";
            throw std::logic_error(ss.str());
        }
    }
    catch(std::exception& e) {
        std::cerr << "Error: " << e.what()
                  << ". Please use -h/--help for correct usage.";
        return 64;  // BSD usage error code.
    }
    catch(...) {
        std::cerr << "Unknown error!" << "\n";
        return 64;  // BSD usage error code.
    }

    /////////////////////////////////////////
    // Load config file and show basic info.
    /////////////////////////////////////////

    auto params = BRKGA::readConfiguration(config_file);
    auto& brkga_params = params.first;
    auto& control_params = params.second;

    // Main algorithm.
    try {
        { // Local scope for start_time.
        auto start_time = system_clock::to_time_t(system_clock::now());

        cout <<
        "------------------------------------------------------\n" <<
        "> Experiment started at " << ctime(&start_time) <<
        "> Instance: " << instance_file <<
        "\n> Configuration: " << config_file <<
        "\n> Algorithm Parameters";
        }

        if(!perform_evolution)
            cout << ">    - Simple multi-start: on (no evolutionary operators)";
        else {
            cout <<
            "\n>   - population_size: " << brkga_params.population_size <<
            "\n>   - elite_percentage: " << brkga_params.elite_percentage <<
            "\n>   - mutants_percentage: " << brkga_params.mutants_percentage <<
            "\n>   - mutants_increase_percentage: " << brkga_params.mutants_increase_percentage <<
            "\n>   - num_elite_parents: " << brkga_params.num_elite_parents <<
            "\n>   - total_parents: " << brkga_params.total_parents <<
            "\n>   - bias_type: " << brkga_params.bias_type <<
            "\n>   - num_independent_populations: " <<
                brkga_params.num_independent_populations <<
            "\n>   - pr_number_pairs: " << brkga_params.pr_number_pairs <<
            "\n>   - pr_minimum_distance: " <<
                brkga_params.pr_minimum_distance <<
            "\n>   - pr_type: " << brkga_params.pr_type <<
            "\n>   - pr_selection: " << brkga_params.pr_selection <<
            "\n>   - alpha_block_size: " << brkga_params.alpha_block_size <<
            "\n>   - pr_percentage: " << brkga_params.pr_percentage <<
            "\n>   - pr_delivery: " << brkga_params.pr_delivery <<
            "\n>   - exchange_interval: " << control_params.exchange_interval <<
            "\n>   - num_exchange_indivuduals: " <<
                control_params.num_exchange_indivuduals <<
            "\n>   - reset_interval: " << control_params.reset_interval;
        }

        cout <<
        "\n> Seed: " << seed <<
        "\n> Stop rule: " <<
            (stop_rule == StopRule::GENERATIONS ? "Generations" :
             (stop_rule == StopRule::TARGET ? "Target" : "Improvement")) <<
        "\n> Stop argument: " << stop_arg <<
        "\n> Maximum time (s): " << max_time <<
        "\n> Number of parallel threads for decoding: " << num_threads <<
        "\n------------------------------------------------------" << endl;


        ////////////////////////////////////////
        // Load instance and adjust BRKGA parameters
        ////////////////////////////////////////
        log("Reading VRPODTW data...");

        auto instance = VRPODTW_Instance(instance_file);
        cout << "Number of nodes: " << instance.num_nodes << endl;

        ////////////////////////////////////////
        // Build the BRKGA data structures
        ////////////////////////////////////////
        log("Building BRKGA...");

        brkga_params.population_size = brkga_params.population_size * instance.num_nodes;
        cout << "New population size: " << brkga_params.population_size << endl;

        // We may choose the decoder
        typedef VRPODTW_Decoder LocalDecoder;
        VRPODTW_Decoder decoder(instance, brkga_params.pr_delivery);

        BRKGA::BRKGA_MP_IPR<LocalDecoder> algorithm(
                decoder, BRKGA::Sense::MINIMIZE, seed, instance.num_nodes,
                brkga_params, num_threads, perform_evolution);


        ////////////////////////////////////////
        // Initialize BRKGA
        ////////////////////////////////////////

        // NOTE: don't forget to initialize the algorithm.
        log("Initializing BRKGA...");
        algorithm.initialize();

        ////////////////////////////////////////
        // Evolving
        ////////////////////////////////////////

        // Distance functor for path-relink.
        std::shared_ptr<BRKGA::DistanceFunctionBase> dist_func;
        if(brkga_params.pr_type == BRKGA::PathRelinking::Type::DIRECT)
            dist_func.reset(new BRKGA::HammingDistance(0.5));
        else
            dist_func.reset(new BRKGA::KendallTauDistance());

        // Optimization info.
        double best_fitness = numeric_limits<double>::max();
        BRKGA::Chromosome best_solution(instance.num_nodes, 0.0);

        unsigned last_update_iteration = 0;
        double last_update_time;
        unsigned large_offset = 0;
        double path_relink_time = 0.0;
        unsigned num_path_relink_calls = 0;
        unsigned restart_wku = 0;
        unsigned num_homogenities = 0;
        unsigned num_best_improvements = 0;
        unsigned num_elite_improvements = 0;
        int cont_mutation = 0;
        int reset_interval = int(population_size * 0.3);
        unsigned num_reset_calls = 0;

        // For control the optimization
        unsigned iteration = 0;
        bool run = true;
        double mutants_pcg = brkga_params.mutants_percentage;

        log("Evolving...");
        cout << "* Iteration | Cost | CurrentTime" << endl;

        const auto start_time = system_clock::now();
        int cont_pr = 0;
        while(run) {
            // Run one iteration
            algorithm.evolve(1, mutants_pcg);

            double fitness = algorithm.getBestFitness();
            if ((best_fitness - fitness) > 0.0) {//data updates for reset and best fitness
                last_update_time = elapsedFrom(start_time);
                best_fitness = fitness;
                restart_wku = num_reset_calls;

                if (iteration - last_update_iteration > M_interval_between_2_improvs)
                    M_interval_between_2_improvs = iteration - last_update_iteration;

                auto update_offset = iteration - last_update_iteration;
                last_update_iteration = iteration;
                if (iteration - last_update_iteration)
                    M_interval_between_2_improvs = iteration;

                if (large_offset < update_offset)
                    large_offset = update_offset;

                best_solution = algorithm.getBestChromosome();

                cout <<
                     "* " << iteration << " | " <<
                     setiosflags(ios::fixed) << setprecision(2) <<
                     best_fitness << " | " <<
                     setiosflags(ios::fixed) << setprecision(2) <<
                     last_update_time <<
                     endl;
            }

            unsigned iters_no_improvement = iteration - last_update_iteration;

            if (iters_no_improvement == reset_interval / 2 &&//START_NEW_MUTANTS
                cont_mutation < 2 &&
                brkga_params.mutants_increase_percentage > 0) {
                cont_mutation += 1;
                mutants_pcg = mutants_pcg + brkga_params.mutants_increase_percentage;
                last_update_iteration = iteration;
            }//END_NEW_MUTANTS

            if (control_params.exchange_interval > 0 && //Start Implicit Path-Relinking (IPR)
                iters_no_improvement > 0 &&
                (iters_no_improvement % control_params.exchange_interval == 0)) {

                num_path_relink_calls += 1;
                cont_pr += 1;

                const auto pr_start_time = system_clock::now();
                auto result = algorithm.pathRelink(
                        brkga_params.pr_type,
                        brkga_params.pr_selection,
                        dist_func,
                        brkga_params.pr_number_pairs,
                        brkga_params.pr_minimum_distance,
                        1, // block_size doesn't not matter for permutation.
                        max_time - elapsedFrom(start_time),
                        brkga_params.pr_percentage,
                        instance.num_nodes - instance.num_drivers);

                const auto pr_time = elapsedFrom(pr_start_time);
                path_relink_time += pr_time;

                using BRKGA::PathRelinking::PathRelinkingResult;
                switch (result) {
                    case PathRelinkingResult::TOO_HOMOGENEOUS:
                        ++num_homogenities;
                        break;

                    case PathRelinkingResult::NO_IMPROVEMENT:
                        break;

                    case PathRelinkingResult::ELITE_IMPROVEMENT:
                        ++num_elite_improvements;
                        break;

                    case PathRelinkingResult::BEST_IMPROVEMENT:
                        ++num_best_improvements;
                        fitness = algorithm.getBestFitness();

                        if ((best_fitness - fitness) > 0.0) {
                            last_update_time = elapsedFrom(start_time);
                            best_fitness = fitness;

                            auto update_offset = iteration -
                                                 last_update_iteration;
                            last_update_iteration = iteration;

                            if (large_offset < update_offset)
                                large_offset = update_offset;

                            best_solution = algorithm.getBestChromosome();

                            cout <<
                                 "* " << iteration << " | " <<
                                 setiosflags(ios::fixed) << setprecision(2) <<
                                 best_fitness << " | " <<
                                 setiosflags(ios::fixed) << setprecision(2) <<
                                 last_update_time <<
                                 endl;
                        }
                        break;

                    default:;
                } // end switch
            } // End if IPR


            if (control_params.reset_interval > 0 &&
                iters_no_improvement > 0 &&//START_NEW_RESET
                (iters_no_improvement % control_params.reset_interval == 0)) {
                mutants_pcg = brkga_params.mutants_percentage;
                cont_mutation = 0;
                num_reset_calls += 1;

                seed = rand() % 27000001;
                algorithm.reset();

                //injectChromosome
                algorithm.injectChromosome(algorithm.getBestChromosome(),
                                           1,
                                           0,
                                           best_fitness);

                cont_pr = 0;
                fitness = algorithm.getBestFitness();

                if ((best_fitness - fitness) > 0.0) {
                    last_update_time = elapsedFrom(start_time);
                    best_fitness = fitness;

                    auto update_offset = iteration -
                                         last_update_iteration;
                    last_update_iteration = iteration;

                    if (large_offset < update_offset)
                        large_offset = update_offset;
                    best_solution = algorithm.getBestChromosome();
                }
            }//END_NEW_RESET

            // Time to stop?
            switch (stop_rule) {
                case StopRule::GENERATIONS:
                    if (iteration == static_cast<unsigned>(stop_arg))
                        run = false;
                    break;

                case StopRule::IMPROVEMENT:
                    if (iters_no_improvement >= static_cast<unsigned>(stop_arg))
                        run = false;
                    break;

                case StopRule::TARGET:
                    if (best_fitness <= stop_arg)
                        run = false;
                    break;
                default:;
            }

            if (elapsedFrom(start_time) > max_time)
                run = false;
            ++iteration;
        }// end while / main loop

        const auto total_num_iterations = iteration;
        const auto total_elapsed_time = elapsedFrom(start_time);

        log("End of optimization");
        cout <<
        "\nTotal number of iterations: " << total_num_iterations <<
        "\nLast update iteration: " << last_update_iteration <<
        "\nTotal optimization time: " << total_elapsed_time <<
        "\nLast update time: " << last_update_time <<
        "\nLarge number of iterations between improvements: " << large_offset <<
        "\nMaximum distance between two improvements: " << M_interval_between_2_improvs <<
        "\nTotal path relink time: " << path_relink_time <<
        "\nTotal path relink calls: " << num_path_relink_calls <<
        "\nTotal reset calls: " << num_reset_calls << "(" <<restart_wku<< ")"
        "\nNumber of homogenities: " << num_homogenities <<
        "\nImprovements in the elite set: " << num_elite_improvements <<
        "\nBest individual improvements: " << num_best_improvements <<
        endl;

        ////////////////////////////////////////
        // Extracting the best solution
        ////////////////////////////////////////

        best_fitness = algorithm.getBestFitness();
        best_solution = algorithm.getBestChromosome();

        cout << setiosflags(ios::fixed) << setprecision(2)
            << best_fitness;

        string instance_name(instance_file);
        size_t pos = instance_name.rfind('/');
        if(pos != string::npos)
                instance_name = instance_name.substr(pos + 1);

        pos = instance_name.rfind('.');
        if(pos != string::npos)
                instance_name = instance_name.substr(0, pos);

        cout <<
        "\n\nInstance,Seed,NumNodes,TotalIterations,TotalTime,"
        "TotalPRTime,PRCalls,NumHomogenities,NumPRImprovElite,"
        "NumPrImprovBest,LargeOffset,LastUpdateIteration,LastUpdateTime,Cost" <<
        endl;

        cout <<
        instance_name << "," <<
        seed << "," <<
        instance.num_nodes << "," <<
        total_num_iterations << "," <<
        setiosflags(ios::fixed) << setprecision(2) <<
        total_elapsed_time << "," <<
        path_relink_time << "," <<
        num_path_relink_calls << "," <<
        num_homogenities << "," <<
        num_elite_improvements << "," <<
        num_best_improvements << "," <<
        large_offset << "," <<
        last_update_iteration << "," <<
        last_update_time << "," <<
        setiosflags(ios::fixed) << setprecision(2) <<
        best_fitness << endl;

        cout.flush();

    //////////////////////////////////////////////

        auto end_time = system_clock::to_time_t(system_clock::now());
        string ss(ctime(&end_time));
        string file_path = "../" + instance_name + "_R" + std::to_string(num_run) + "_k" + std::to_string(control_params.reset_interval)  + ".txt";
        ofstream f(file_path, ios::out);
        f<<
        "> Experiment ended at " << ctime(&end_time) <<
        "> Instance: " << instance_file <<
        "\n> Configuration: " << config_file <<
        "\n> Algorithm Parameters" <<
        "\n>   - population_size: " << brkga_params.population_size <<
        "\n>   - elite_percentage: " << brkga_params.elite_percentage <<
        "\n>   - mutants_percentage: " << brkga_params.mutants_percentage <<
        "\n>   - mutants_increase_percentage: " << brkga_params.mutants_increase_percentage <<
        "\n>   - num_elite_parents: " << brkga_params.num_elite_parents <<
        "\n>   - total_parents: " << brkga_params.total_parents <<
        "\n>   - bias_type: " << brkga_params.bias_type <<
        "\n>   - num_independent_populations: " <<
        brkga_params.num_independent_populations <<
        "\n>   - pr_number_pairs: " << brkga_params.pr_number_pairs <<
        "\n>   - pr_minimum_distance: " <<
        brkga_params.pr_minimum_distance <<
        "\n>   - pr_type: " << brkga_params.pr_type <<
        "\n>   - pr_selection: " << brkga_params.pr_selection <<
        "\n>   - alpha_block_size: " << brkga_params.alpha_block_size <<
        "\n>   - pr_percentage: " << brkga_params.pr_percentage <<
        "\n>   - pr_delivery: " << brkga_params.pr_delivery <<
        "\n>   - exchange_interval: " << control_params.exchange_interval <<
        "\n>   - num_exchange_indivuduals: " <<
        control_params.num_exchange_indivuduals <<
        "\n>   - reset_interval: " << control_params.reset_interval <<

        "\n> Seed: " << seed <<
        "\n> Stop rule: " <<
        (stop_rule == StopRule::GENERATIONS ? "Generations" :
        (stop_rule == StopRule::TARGET ? "Target" : "Improvement")) <<
        "\n> Stop argument: " << stop_arg <<
        "\n> Maximum time (s): " << max_time <<
        "\n> Number of parallel threads for decoding: " << num_threads <<
        "\n------------------------------------------------------" <<
        "\nTotal number of iterations: " << total_num_iterations <<
        "\nLast update iteration: " << last_update_iteration <<
        "\nTotal optimization time: " << total_elapsed_time <<
        "\nLast update time: " << last_update_time <<
        "\nLarge number of iterations between improvements: " << large_offset <<
        "\nMaximum distance between two improvements: " << M_interval_between_2_improvs <<
        "\nTotal path relink time: " << path_relink_time <<
        "\nTotal path relink calls: " << num_path_relink_calls <<
        "\nTotal reset calls: " << num_reset_calls << "(" <<restart_wku<< ")"
        "\nNumber of homogenities: " << num_homogenities <<
        "\nImprovements in the elite set: " << num_elite_improvements <<
        "\nBest individual improvements: " << num_best_improvements <<

        "\n% Best "
        << setiosflags(ios::fixed) << setprecision(2)
        << best_fitness;

        f<<"\nInstance,Seed,NumNodes,TotalIterations,TotalTime,"
        "TotalPRTime,PRCalls,NumHomogenities,NumPRImprovElite,"
        "NumPrImprovBest,LargeOffset,LastUpdateIteration,LastUpdateTime,Cost" <<
        "\n" << instance_name << "," <<
        seed << "," <<
        instance.num_nodes << "," <<
        total_num_iterations << "," <<
        setiosflags(ios::fixed) << setprecision(2) <<
        total_elapsed_time << "," <<
        path_relink_time << "," <<
        num_path_relink_calls << "," <<
        num_homogenities << "," <<
        num_elite_improvements << "," <<
        num_best_improvements << "," <<
        large_offset << "," <<
        last_update_iteration << "," <<
        last_update_time << "," <<
        setiosflags(ios::fixed) << setprecision(2) <<
        best_fitness << endl;
        f.close();
    //////////////////////////////////////////////
    // Writing the output to a txt file - END
    //////////////////////////////////////////////
    }
    catch(exception& e) {
        cerr << "\n***********************************************************"
             << "\n****  Exception Occured: " << e.what()
             << "\n***********************************************************"
             << endl;
        return 70; // BSD software internal error code
    }
    return 0;
}
