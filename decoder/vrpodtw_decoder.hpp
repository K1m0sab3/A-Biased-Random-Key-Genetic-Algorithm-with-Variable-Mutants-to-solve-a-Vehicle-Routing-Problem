/******************************************************************************
 * vrpodtw_decoder.hpp: interface for VRPODTW_Decoder class.
 *****************************************************************************/

#ifndef VRPODTW_DECODER_HPP_
#define VRPODTW_DECODER_HPP_
#define COMPENSATION 1.1


#include "../reading_instances/vrpodtw_instance.hpp"
#include "../../../brkga_mp_ipr/chromosome.hpp"

class VRPODTW_Decoder {
public:
    /** \brief Default Constructor.
     * \param instance VRPODTW instance.
     */
    VRPODTW_Decoder(const VRPODTW_Instance& instance, int pr_delivery_value);


/** \brief Given a chromossome, builds a set of paths.
     *
     * \param chromosome A vector of doubles represent a problem solution.
     * \param rewrite Indicates if the chromosome must be rewritten. Not used
     *                this decoder, but keep due to API requirements.
     * \return the cost of the paths.
     */
double decode(BRKGA::Chromosome& chromosome, bool rewrite, double currentBestFitness);
public:
    /// A reference to a VRPODTW instance.
    const VRPODTW_Instance& instance;
    std::vector<std::vector<int>> paths;
    int pr_delivery_value;

};

#endif // VRPODTW_DECODER_HPP_
