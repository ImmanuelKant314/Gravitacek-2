#include <stdexcept>
#include <cmath>
#include <iostream>

#include "gravitacek2/mymath.hpp"
#include "gravitacek2/geomotion/spacetimes.hpp"

namespace gr2
{
    ReissnerNordstromMPW::ReissnerNordstromMPW(const real& M)
    {
        this->M = M;
    };

    ReissnerNordstromMPW::~ReissnerNordstromMPW()
    {

    };

    void ReissnerNordstromMPW::calculate_N(const real* y) 
    {
        // coordinates
        real rho = y[RHO];
        real z = y[Z];
        real d = sqrtl(rho*rho + z*z);
        real d_inv = 1.0/d;

        // lapse
        this->N_inv = 1.0 + M*d_inv;
        this->N = 1.0/N_inv;
    };

    void ReissnerNordstromMPW::calculate_N1(const real* y)
    {
        // coordinates
        real rho = y[RHO];
        real z = y[Z];
        real d = sqrtl(rho*rho + z*z);
        real d_inv = 1.0/d;
        real M_d_inv = 1.0/(M+d);

        // lapse
        this->N_inv = 1.0 + M*d_inv;
        this->N = 1.0/N_inv;

        // derivative of lapse
        this->N_rho = M*rho*M_d_inv*M_d_inv*d_inv;
        this->N_z = M*z*M_d_inv*M_d_inv*d_inv;
    };

    void ReissnerNordstromMPW::calculate_N2(const real* y)
    {
        // coordinates
        real rho = y[RHO];
        real z = y[Z];

        // Second derivatives of potential - functions
        auto N_rho_func = [&z, this](real rho)
        {
            real y[] = {0, 0, 0, 0};
            y[RHO] = rho;
            y[Z] = z;
            this->calculate_N1(y);
            return this->get_N_rho();
        };

        auto N_z_func = [&rho, this](real z)
        {
            real y[] = {0, 0, 0, 0};
            y[RHO] = rho;
            y[Z] = z;
            this->calculate_N1(y);
            return this->get_N_z();
        };

        auto N_z_func_rho = [&z, this](real rho)
        {
            real y[] = {0, 0, 0, 0};
            y[RHO] = rho;
            y[Z] = z;
            this->calculate_N1(y);
            return this->get_N_z();
        };

        // Second derivatives of lapse
        this->N_rhorho = gr2::richder<5>(N_rho_func, rho, 0.1, 1e-10);
        this->N_zz = gr2::richder<5>(N_z_func, z, 0.1, 1e-10);
        this->N_rhoz = gr2::richder<5>(N_z_func_rho, rho, 0.1, 1e-10);

        // Prepare variables
        real d = sqrtl(rho*rho + z*z);
        real d_inv = 1.0/d;
        real Md = (M+d);
        real Md_inv = 1.0/(M+d);
        real Md_inv3 = Md_inv*Md_inv*Md_inv;

        // lapse
        this->N_inv = 1.0 + M*d_inv;
        this->N = 1.0/N_inv;

        // derivative of lapse
        this->N_rho = M*rho*Md_inv*Md_inv*d_inv;
        this->N_z = M*z*Md_inv*Md_inv*d_inv;
    };
}