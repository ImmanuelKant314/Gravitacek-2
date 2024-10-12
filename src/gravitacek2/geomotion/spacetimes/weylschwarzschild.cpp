#include <cmath>

#include "gravitacek2/geomotion/spacetimes.hpp"

namespace gr2
{
    WeylSchwarzschild::WeylSchwarzschild(real M, bool lambda_exact) : Weyl(lambda_exact)
    {
        this->M = M;
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

    void WeylSchwarzschild::calculate_lambda(const real *y)
    {
        real rho = y[RHO];
        real z = y[Z];
        real d1 = sqrtl(rho*rho + (z-M)*(z-M));
        real d2 = sqrtl(rho*rho + (z+M)*(z+M));

        real Sigma = d1*d2;

        this->lambda = 0.5*logl(((d1+d2)*(d1+d2)-4*M*M)/(4*Sigma));
    }
}