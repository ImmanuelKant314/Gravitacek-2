#include <stdexcept>

#include "gravitacek2/geomotion/majumadpapapetrouweyl.hpp"

namespace gr2
{
    CombinedMPW::CombinedMPW(std::vector<std::shared_ptr<MajumdarPapapetrouWeyl>> sources): sources(sources)
    {

    };

    CombinedMPW::~CombinedMPW()
    {

    };

    void CombinedMPW::calculate_N_inv(const real* y)
    {
        this->N_inv = 1;
        for (auto &s : sources)
        {
            s->calculate_N_inv(y);
            N_inv += s->get_N_inv() - 1;
        }
    };

    void CombinedMPW::calculate_N_inv1(const real* y)
    {
        this->N_inv = 1;
        this->N_inv_rho = 0;
        this->N_inv_z = 0;
        for (auto &s : sources)
        {
            s->calculate_N_inv1(y);
            N_inv += s->get_N_inv() - 1;
            N_inv_rho += s->get_N_inv_rho();
            N_inv_z += s->get_N_inv_z();
        }
    };

    void CombinedMPW::calculate_N_inv2(const real* y)
    {
        this->N_inv = 1;
        this->N_inv_rho = 0;
        this->N_inv_z = 0;
        this->N_inv_rhorho = 0;
        this->N_inv_rhoz = 0;
        this->N_inv_zz = 0;

        for (auto &s : sources)
        {
            s->calculate_N_inv2(y);
            N_inv += s->get_N_inv() - 1;
            N_inv_rho += s->get_N_inv_rho();
            N_inv_z += s->get_N_inv_z();
            N_inv_rhorho += s->get_N_inv_rhorho();
            N_inv_rhoz += s->get_N_inv_rhoz();
            N_inv_zz += s->get_N_inv_zz();
        }
    }

}