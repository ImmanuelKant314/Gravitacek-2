#include <cmath>
#include <stdexcept>

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
        real rho = y[RHO];
        real z = y[Z];
        real d1 = sqrtl(rho*rho + (z-M)*(z-M));
        real d2 = sqrtl(rho*rho + (z+M)*(z+M));

        real r = 0.5 * (d1 + d2) + M;

        this->nu = 0.5*logl((d1+d2-2*M)/(d1+d2+2*M));
        this->nu_rho = M*rho/(2*r*(r-2*M))*(1.0/d1 + 1.0/d2);
        this->nu_z = M/(2*r*(r-2*M))*((z-M)/d1 + (z+M)/d2);
    }
}