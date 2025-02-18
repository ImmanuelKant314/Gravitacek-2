#include <cmath>
#include <stdexcept>

#include "gravitacek2/mymath.hpp"
#include "gravitacek2/geomotion/spacetimes.hpp"

namespace gr2
{
    void WeylSchwarzschild::calculate_lambda_exact(const real *y)
    {
        real rho = y[RHO];
        real z = y[Z];
        real d1 = sqrtl(rho*rho + (z-M)*(z-M));
        real d2 = sqrtl(rho*rho + (z+M)*(z+M));

        real Sigma = d1*d2;

        this->lambda = 0.5*logl(((d1+d2)*(d1+d2)-4*M*M)/(4*Sigma));
    }

    void WeylSchwarzschild::calculate_lambda_integral(const real *y)
    {
        calculate_lambda_from_inf_to_z(y, 1e-15);
    }

    WeylSchwarzschild::WeylSchwarzschild(real M, LambdaEvaluation init, LambdaEvaluation run) : Weyl(init, run)
    {
        this->M = M;
    }

    void WeylSchwarzschild::calculate_lambda_init(real const* y)
    {
        switch (this->lambda_eval_init)
        {
        case LambdaEvaluation::exact:
            this->calculate_lambda_exact(y);
            break;
        case LambdaEvaluation::integral:
            this->calculate_lambda_integral(y);
            break;
        default:
            throw std::runtime_error("Calculating lambda this way is not possible");
            break;
        }
    }

    void WeylSchwarzschild::calculate_lambda_run(real const* y)
    {
        switch (lambda_eval_run)
        {
        case LambdaEvaluation::exact:
            this->calculate_lambda_exact(y);
            break;
        case LambdaEvaluation::diff:
            this->calculate_lambda_diff(y);
            break;
        default:
            throw std::runtime_error("Calculating lambda this way is not possible");
            break;
        }
    }

    void WeylSchwarzschild::calculate_nu(const real* y) 
    {
        real rho = y[RHO];
        real z = y[Z];
        real d1 = sqrtl(rho*rho + (z-M)*(z-M));
        real d2 = sqrtl(rho*rho + (z+M)*(z+M));

        this->nu = 0.5*logl((d1+d2-2*M)/(d1+d2+2*M));
    }

    void WeylSchwarzschild::calculate_nu1(const real* y)
    {
        real rho = y[RHO];
        real z = y[Z];
        real d1 = sqrtl(rho*rho + (z-M)*(z-M));
        real d2 = sqrtl(rho*rho + (z+M)*(z+M));

        real r = 0.5 * (d1 + d2) + M;

        this->nu = 0.5*logl((d1+d2-2*M)/(d1+d2+2*M));
        this->nu_rho = M*rho/(2*r*(r-2*M))*(1.0/d1 + 1.0/d2);
        this->nu_z = M/(2*r*(r-2*M))*((z-M)/d1 + (z+M)/d2);
    }

    void WeylSchwarzschild::calculate_nu2(const real* y)
    {
        // Coordinates
        real rho = y[RHO];
        real z = y[Z];

        // Second derivatives of potential
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

        // Prepare calculations
        real d1 = sqrtl(rho*rho + (z-M)*(z-M));
        real d2 = sqrtl(rho*rho + (z+M)*(z+M));
        real r = 0.5 * (d1 + d2) + M;

        // Potential
        this->nu = 0.5*logl((d1+d2-2*M)/(d1+d2+2*M));

        // First derivatives of potential
        this->nu_rho = M*rho/(2*r*(r-2*M))*(1.0/d1 + 1.0/d2);
        this->nu_z = M/(2*r*(r-2*M))*((z-M)/d1 + (z+M)/d2);
    }
}