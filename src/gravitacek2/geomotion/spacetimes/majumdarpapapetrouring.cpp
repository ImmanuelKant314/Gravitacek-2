#include <stdexcept>
#include <cmath>

#include "gravitacek2/mymath.hpp"
#include "gravitacek2/geomotion/spacetimes.hpp"

namespace gr2
{
    MajumdarPapapetrouRing::MajumdarPapapetrouRing(const real &M, const real &b):M(M), b(b)
    {

    };

    MajumdarPapapetrouRing::~MajumdarPapapetrouRing()
    {

    };

    void MajumdarPapapetrouRing::calculate_N_inv(const real* y)
    {
        real rho = y[RHO], z = y[Z];
        real K, E;
        real l1 = sqrtl((rho-b)*(rho-b) + z*z);
        real l2 = sqrtl((rho+b)*(rho+b) + z*z);
        real k = sqrtl(4*rho*b/(l2*l2));
        elliptic_KE(k, K, E);

        this->N_inv = 1 + 2*M*K/(pi*l2);
    };


    void MajumdarPapapetrouRing::calculate_N_inv1(const real* y)
    {
        real rho = y[RHO], z = y[Z];
        real K, E;
        real l1 = sqrtl((rho-b)*(rho-b) + z*z);
        real l2 = sqrtl((rho+b)*(rho+b) + z*z);
        real k = sqrtl(4*rho*b/(l2*l2));
        elliptic_KE(k, K, E);

        this->N_inv = 1 + 2*M*K/(pi*l2);
        this->N_inv_rho = -M*(l1*l1*K-(b*b + z*z - rho*rho)*E)/(pi*rho*l1*l1*l2);
        this->N_inv_z = -2*M*E*z/(pi*l1*l1*l2);
    };


    void MajumdarPapapetrouRing::calculate_N_inv2(const real* y)
    {
        real rho = y[RHO], z = y[Z];

        auto N_inv_rho_func = [&z, this](real rho)
        {
            real y[] = {0, 0, 0, 0};
            y[RHO] = rho;
            y[Z] = z;
            this->calculate_N_inv1(y);
            return this->get_N_inv_rho();
        };
        
        auto N_inv_z_func = [&rho, this](real z)
        {
            real y[] = {0, 0, 0, 0};
            y[RHO] = rho;
            y[Z] = z;
            this->calculate_N_inv1(y);
            return this->get_N_inv_z();
        };

        auto N_inv_z_func_rho = [&z, this](real rho)
        {
            real y[] = {0, 0, 0, 0};
            y[RHO] = rho;
            y[Z] = z;
            this->calculate_N_inv1(y);
            return this->get_N_inv_z();
        };

        // Second derivatives of N_inv
        this->N_inv_rhorho = gr2::richder<5>(N_inv_rho_func, rho, 0.1, 1e-10);
        this->N_inv_zz = gr2::richder<5>(N_inv_z_func, z, 0.1, 1e-10);
        this->N_inv_rhoz = gr2::richder<5>(N_inv_z_func_rho, rho, 0.1, 1e-10);

        real K, E;
        real l1 = sqrtl((rho-b)*(rho-b) + z*z);
        real l2 = sqrtl((rho+b)*(rho+b) + z*z);
        real k = sqrtl(4*rho*b/(l2*l2));
        elliptic_KE(k, K, E);

        this->N_inv = 1 + 2*M*K/(pi*l2);
        this->N_inv_rho = -M*(l1*l1*K-(b*b + z*z - rho*rho)*E)/(pi*rho*l1*l1*l2);
        this->N_inv_z = -2*M*E*z/(pi*l1*l1*l2);
    };
}