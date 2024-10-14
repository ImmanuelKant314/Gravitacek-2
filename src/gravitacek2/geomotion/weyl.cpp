#include <cmath>
#include <stdexcept>

#include "gravitacek2/geomotion/weyl.hpp"
#include "gravitacek2/mymath.hpp"

namespace gr2
{
    void Weyl::calculate_lambda_from_inf_to_z(const real* y, const real& eps)
    {
        // add variables
        real a = 0.5;
        real rho = y[RHO];
        real z = y[Z];
        real b = 1.0/(z+1); 
        const int I_MAX = 20;

        // create integrated function
        auto integrated_function = [&rho, this](real x)
        {
            real y[] = {0, 0, 0, 0};
            y[RHO] = rho;
            y[Z] = 1.0/x-1;
            this->calculate_nu1(y);
            return -(2*rho*nu_rho*nu_z)/(x*x);
        };

        // find value of a
        int i;
        for (i = 0; i < I_MAX; i++)
        {
            real val = integrated_function(a);
            if (fabsl(val)*a < eps)
                break;
            a *= 0.5;
        }

        if (i == I_MAX)
            throw std::runtime_error("in calculating lambda lower bound was not found");

        // calculate value lambda
        this->lambda = romb<5>(integrated_function, a, b, eps);
    }

    Weyl::Weyl(LambdaEvaluation init, LambdaEvaluation run) : GeoMotion(4, run==LambdaEvaluation::diff?9:8)
    {
        // ways of calculating lambda
        this->lambda_eval_init = init;
        this->lambda_eval_run = run;

        // index of lambda
        this->lambda_index = LAMBDA;
    }

    Weyl::~Weyl()
    {

    }

    real Weyl::get_nu() const
    {
        return this->nu;
    }

    real Weyl::get_nu_rho() const
    {
        return this->nu_rho;
    }

    real Weyl::get_nu_z() const
    {
        return this->nu_z;
    }

    real Weyl::get_nu_rhorho() const
    {
        return this->nu_rhorho;
    }

    real Weyl::get_nu_rhoz() const
    {
        return this->nu_rhoz;
    }

    real Weyl::get_nu_zz() const
    {
        return this->nu_zz;
    }

    void Weyl::calculate_lambda_diff(const real* y)
    {
        this->lambda = y[lambda_index];
    }

    real Weyl::get_lambda() const
    {
        return this->lambda;
    }

    void Weyl::set_lambda_index(const int& lambda_index)
    {
        this->lambda_index = lambda_index;
    }

    int Weyl::get_lambda_index() const
    {
        return this->lambda_index;
    }

    void Weyl::calculate_metric(const real *y)
    {
        if(!necessary_calculate(y, y_m, dim))
            return;

        real t = y[T];
        real phi = y[PHI];
        real rho = y[RHO];
        real z = y[Z];

        this->calculate_lambda_run(y); 
        this->calculate_nu(y);

        real exp_2nu = expl(2*nu);
        real exp_2nu_inv = expl(-2*nu);
        real exp_2lambda_2nu_inv = expl(2*lambda-2*nu);

        metric[T][T] = -exp_2nu;
        metric[PHI][PHI] = rho*rho*exp_2nu_inv;
        metric[RHO][RHO] = exp_2lambda_2nu_inv;
        metric[Z][Z] = exp_2lambda_2nu_inv;
    };

    void Weyl::calculate_christoffel_symbols(const real *y)
    {
        if(!necessary_calculate(y, y_c, dim))
            return;

        real t = y[T];
        real phi = y[PHI];
        real rho = y[RHO];
        real z = y[Z];

        this->calculate_lambda_run(y);
        this->calculate_nu1(y);

        real exp_4nu = expl(4*nu);
        real exp_2lambda_inv = expl(-2*lambda);
        this->lambda_rho = rho*(nu_rho*nu_rho - nu_z*nu_z);
        this->lambda_z = 2*rho*nu_rho*nu_z;

        christoffel_symbols[T][T][RHO] = nu_rho;
        christoffel_symbols[T][RHO][T] = christoffel_symbols[T][T][RHO];
        christoffel_symbols[T][T][Z] = nu_z;
        christoffel_symbols[T][Z][T] = christoffel_symbols[T][T][Z];
        christoffel_symbols[PHI][PHI][RHO] = -nu_rho + 1.0/rho;
        christoffel_symbols[PHI][RHO][PHI] = christoffel_symbols[PHI][PHI][RHO];
        christoffel_symbols[PHI][PHI][Z] = -nu_z;
        christoffel_symbols[PHI][Z][PHI] = christoffel_symbols[PHI][PHI][Z];
        christoffel_symbols[RHO][T][T] = exp_2lambda_inv*exp_4nu*nu_rho;
        christoffel_symbols[RHO][PHI][PHI] = rho*(rho*nu_rho-1)*exp_2lambda_inv;
        christoffel_symbols[RHO][RHO][RHO] = lambda_rho-nu_rho;
        christoffel_symbols[RHO][RHO][Z] = lambda_z - nu_z;
        christoffel_symbols[RHO][Z][RHO] = christoffel_symbols[RHO][RHO][Z];
        christoffel_symbols[RHO][Z][Z] = -lambda_rho+nu_rho;
        christoffel_symbols[Z][T][T] = exp_2lambda_inv*exp_4nu*nu_z;
        christoffel_symbols[Z][PHI][PHI] = rho*rho*exp_2lambda_inv*nu_z;
        christoffel_symbols[Z][RHO][RHO] = -lambda_z + nu_z;
        christoffel_symbols[Z][RHO][Z] = lambda_rho - nu_rho;
        christoffel_symbols[Z][Z][RHO] = christoffel_symbols[Z][RHO][Z];
        christoffel_symbols[Z][Z][Z] = lambda_z - nu_z;
    };

    void Weyl::calculate_riemann_tensor(const real *y)
    {

    };
}