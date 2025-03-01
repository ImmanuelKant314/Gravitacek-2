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
        if(!necessary_calculate(y, y_r, dim))
            return;

        real t = y[T];
        real phi = y[PHI];
        real rho = y[RHO];
        real z = y[Z];

        this->calculate_lambda_run(y);
        this->calculate_nu2(y);

        real exp_4nu = expl(4*nu);
        real exp_2lambda_inv = expl(-2*lambda);
        real rho_inv = 1.0/rho;

        // First derivatives of lambda
        this->lambda_rho = rho*(nu_rho*nu_rho - nu_z*nu_z);
        this->lambda_z = 2*rho*nu_rho*nu_z;

        // Second derivatives of lambda
        this->lambda_rhorho = (nu_rho*nu_rho - nu_z*nu_z) + 2*rho*(nu_rho*nu_rhorho - nu_z*nu_rhoz);
        this->lambda_rhoz = 2*rho*(nu_rho*nu_rhoz - nu_z*nu_zz);
        this->lambda_zz = 2*rho*(nu_rhoz*nu_z + nu_rho*nu_zz);

        // Riemann tensor
        riemann_tensor[T][PHI][T][PHI] = rho*(rho*nu_z*nu_z + (rho*nu_rho-1)*nu_rho)*exp_2lambda_inv;
        riemann_tensor[T][PHI][PHI][T] = -riemann_tensor[T][PHI][T][PHI];
        riemann_tensor[T][RHO][T][RHO] = lambda_rho*nu_rho - lambda_z*nu_z - 2*nu_rho*nu_rho - nu_rhorho + nu_z*nu_z;
        riemann_tensor[T][RHO][RHO][T] = -riemann_tensor[T][RHO][T][RHO];
        riemann_tensor[T][RHO][T][Z] = lambda_rho*nu_z + lambda_z*nu_rho - 3*nu_rho*nu_z - nu_rhoz;
        riemann_tensor[T][RHO][Z][T] = -riemann_tensor[T][RHO][T][Z];
        riemann_tensor[T][Z][T][RHO] = riemann_tensor[T][RHO][T][Z];
        riemann_tensor[T][Z][RHO][T] = -riemann_tensor[T][Z][T][RHO];
        riemann_tensor[T][Z][T][Z] = -lambda_rho*nu_rho + lambda_z*nu_z + nu_rho*nu_rho - 2*nu_z*nu_z - nu_zz;
        riemann_tensor[T][Z][Z][T] = -riemann_tensor[T][Z][T][Z];
        riemann_tensor[PHI][T][T][PHI] = rho_inv*(rho*nu_rho*nu_rho + rho*nu_z*nu_z - nu_rho)*exp_4nu*exp_2lambda_inv;
        riemann_tensor[PHI][T][PHI][T] = -riemann_tensor[PHI][T][T][PHI];
        riemann_tensor[PHI][RHO][PHI][RHO] = -lambda_rho*nu_rho + lambda_z*nu_z + nu_rhorho - nu_z*nu_z + rho_inv*(lambda_rho + nu_rho);
        riemann_tensor[PHI][RHO][RHO][PHI] = -riemann_tensor[PHI][RHO][PHI][RHO];
        riemann_tensor[PHI][RHO][PHI][Z] = -lambda_rho*nu_z - lambda_z*nu_rho + nu_rho*nu_z + nu_rhoz + rho_inv*lambda_z;
        riemann_tensor[PHI][RHO][Z][PHI] = -riemann_tensor[PHI][RHO][PHI][Z]; 
        riemann_tensor[PHI][Z][PHI][RHO] = riemann_tensor[PHI][RHO][PHI][Z];
        riemann_tensor[PHI][Z][RHO][PHI] = -riemann_tensor[PHI][Z][PHI][RHO];
        riemann_tensor[PHI][Z][PHI][Z] = lambda_rho*nu_rho-lambda_z*nu_z - nu_rho*nu_rho + nu_zz + rho_inv*(-lambda_rho + nu_rho);
        riemann_tensor[PHI][Z][Z][PHI] = -riemann_tensor[PHI][Z][PHI][Z];
        riemann_tensor[RHO][T][T][RHO] = (lambda_rho*nu_rho - lambda_z*nu_z - 2*nu_rho*nu_rho - nu_rhorho + nu_z*nu_z)*exp_4nu*exp_2lambda_inv;
        riemann_tensor[RHO][T][RHO][T] = -riemann_tensor[RHO][T][T][RHO];
        riemann_tensor[RHO][T][T][Z] = (lambda_rho*nu_z + lambda_z*nu_rho-3*nu_rho*nu_z-nu_rhoz)*exp_4nu*exp_2lambda_inv;
        riemann_tensor[RHO][T][Z][T] = -riemann_tensor[RHO][T][T][Z];
        riemann_tensor[RHO][PHI][PHI][RHO] = rho*(rho*(lambda_rho*nu_rho-lambda_z*nu_z-nu_rhorho+nu_z*nu_z)-lambda_rho-nu_rho)*exp_2lambda_inv;
        riemann_tensor[RHO][PHI][RHO][PHI] = -riemann_tensor[RHO][PHI][PHI][RHO];
        riemann_tensor[RHO][PHI][PHI][Z] = rho*(rho*(lambda_rho*nu_z + lambda_z*nu_rho-nu_rho*nu_z - nu_rhoz) - lambda_z)*exp_2lambda_inv;
        riemann_tensor[RHO][PHI][Z][PHI] = -riemann_tensor[RHO][PHI][PHI][Z];
        riemann_tensor[RHO][Z][RHO][Z] = -lambda_rhorho - lambda_zz + nu_rhorho + nu_zz;
        riemann_tensor[RHO][Z][Z][RHO] = -riemann_tensor[RHO][Z][RHO][Z];
        riemann_tensor[Z][T][T][RHO] = (lambda_rho*nu_z + lambda_z*nu_rho - 3*nu_rho*nu_z - nu_rhoz)*exp_4nu*exp_2lambda_inv;
        riemann_tensor[Z][T][RHO][T] = -riemann_tensor[Z][T][T][RHO];
        riemann_tensor[Z][T][T][Z] = (-lambda_rho*nu_rho + lambda_z*nu_z + nu_rho*nu_rho - 2*nu_z*nu_z - nu_zz)*exp_4nu*exp_2lambda_inv;
        riemann_tensor[Z][T][Z][T] = -riemann_tensor[Z][T][T][Z];
        riemann_tensor[Z][PHI][PHI][RHO] = rho*(rho*(lambda_rho*nu_z + lambda_z*nu_rho-nu_rho*nu_z - nu_rhoz) - lambda_z)*exp_2lambda_inv;
        riemann_tensor[Z][PHI][RHO][PHI] = -riemann_tensor[Z][PHI][PHI][RHO];
        riemann_tensor[Z][PHI][PHI][Z] = rho*(rho*(-lambda_rho*nu_rho+lambda_z*nu_z + nu_rho*nu_rho - nu_zz) + lambda_rho - nu_rho)*exp_2lambda_inv;
        riemann_tensor[Z][PHI][Z][PHI] = -riemann_tensor[Z][PHI][PHI][Z];
        riemann_tensor[Z][RHO][RHO][Z] = lambda_rhorho + lambda_zz - nu_rhorho - nu_zz;
        riemann_tensor[Z][RHO][Z][RHO] = -riemann_tensor[Z][RHO][RHO][Z];
    };


    void Weyl::function(const real &t, const real y[], real dydt[])
    {
        this->GeoMotion::function(t, y, dydt);
        if (this->lambda_eval_run == gr2::diff)
            dydt[this->lambda_index] = this->lambda_rho*y[URHO] + this->lambda_z*y[UZ];
    }
}