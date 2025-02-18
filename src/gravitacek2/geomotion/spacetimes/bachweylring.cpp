#include <cmath>
#include <stdexcept>

#include "gravitacek2/geomotion/spacetimes.hpp"
#include "gravitacek2/mymath.hpp"

namespace gr2
{
    void BachWeylRing::calculate_lambda_integral(const real *y)
    {
        calculate_lambda_from_inf_to_z(y, 1e-15);
    }

    BachWeylRing::BachWeylRing(real M, real b, LambdaEvaluation init, LambdaEvaluation run) : Weyl(init, run)
    {
        this->M = M;
        this->b = b;
    }

    void BachWeylRing::calculate_lambda_init(real const* y)
    {
        switch (this->lambda_eval_init)
        {
        case LambdaEvaluation::integral:
            this->calculate_lambda_integral(y);
            break;
        default:
            throw std::runtime_error("Calculating lambda this way is not possible");
            break;
        }
    }

    void BachWeylRing::calculate_lambda_run(real const* y)
    {
        switch (this->lambda_eval_run)
        {
        case LambdaEvaluation::diff:
            this->calculate_lambda_diff(y);
            break;
        default:
            throw std::runtime_error("Calculating lambda this way is not possible");
            break;
        }
    }

    void BachWeylRing::calculate_nu(const real* y)
    {
        real rho = y[RHO], z = y[Z];
        real K, E;
        real l1 = sqrtl((rho-b)*(rho-b) + z*z);
        real l2 = sqrtl((rho+b)*(rho+b) + z*z);
        real k = sqrtl(4*rho*b/(l2*l2));
        elliptic_KE(k, K, E);
        
        this->nu = -2*M*K/(pi*l2);
    }

    void BachWeylRing::calculate_nu1(const real* y)
    {
        real rho = y[RHO], z = y[Z];
        real K, E;
        real l1 = sqrtl((rho-b)*(rho-b) + z*z);
        real l2 = sqrtl((rho+b)*(rho+b) + z*z);
        real k = sqrtl(4*rho*b/(l2*l2));
        elliptic_KE(k, K, E);
        
        this->nu = -2*M*K/(pi*l2);
        this->nu_rho = M*(l1*l1*K-(b*b + z*z - rho*rho)*E)/(pi*rho*l1*l1*l2);
        this->nu_z = 2*M*E*z/(pi*l1*l1*l2);
    }

    void BachWeylRing::calculate_nu2(const real* y)
    {
        // Coordinates 
        real rho = y[RHO], z = y[Z];

        // Second derivatives of potential - functions
        auto nu_rho_func = [&z, this](real rho)
        {
            real y[] = {0, 0, 0, 0};
            y[RHO] = rho;
            y[Z] = z;
            this->calculate_nu1(y);
            return this->get_nu_rho();
        };

        auto nu_z_func = [&rho, this](real z)
        {
            real y[] = {0, 0, 0, 0};
            y[RHO] = rho;
            y[Z] = z;
            this->calculate_nu1(y);
            return this->get_nu_z();
        };

        auto nu_z_func_rho = [&z, this](real rho)
        {
            real y[] = {0, 0, 0, 0};
            y[RHO] = rho;
            y[Z] = z;
            this->calculate_nu1(y);
            return this->get_nu_z();
        };

        // Second derivatives of potential
        this->nu_rhorho = gr2::richder<5>(nu_rho_func, rho, 0.1, 1e-10);
        this->nu_zz = gr2::richder<5>(nu_z_func, z, 0.1, 1e-10);
        this->nu_rhoz = gr2::richder<5>(nu_z_func_rho, rho, 0.1, 1e-10);

        // Prepare calculation
        real K, E;
        real l1 = sqrtl((rho-b)*(rho-b) + z*z);
        real l2 = sqrtl((rho+b)*(rho+b) + z*z);
        real k = sqrtl(4*rho*b/(l2*l2));
        elliptic_KE(k, K, E);
        
        // Potential
        this->nu = -2*M*K/(pi*l2);

        // First derivatives of potential
        this->nu_rho = M*(l1*l1*K-(b*b + z*z - rho*rho)*E)/(pi*rho*l1*l1*l2);
        this->nu_z = 2*M*E*z/(pi*l1*l1*l2);
    }

}