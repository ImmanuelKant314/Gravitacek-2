#include <cmath>
#include <stdexcept>

#include "gravitacek2/geomotion/spacetimes.hpp"
#include "gravitacek2/mymath.hpp"

namespace gr2
{
    void InvertedKuzminToomreDisk::calculate_lambda_integral(const real *y)
    {
        calculate_lambda_from_inf_to_z(y, 1e-15);
    }

    InvertedKuzminToomreDisk::InvertedKuzminToomreDisk(int n, real M, real b, LambdaEvaluation init, LambdaEvaluation run): Weyl(init, run)
    {
        // save values
        this->n = n;
        this->M = M;
        this->b = b;

        // prepare arrays
        this->B = nullptr;
        this->P0 = nullptr;
        this->P1 = nullptr;

        // create arrays
        this->B = new real[n+1];
        this->P0 = new real[n+1];
        this->P1 = new real[n+1];

        // calculate normalization
        real factorial = 1;
        real double_fact = 1;
        real factorial_half = 1;
        for (int i = 1; i <= n; i++)
        {
            double_fact *= 2*i+1;
            factorial *= i;
            factorial_half *= (real)i+0.5;
        }
        N = -factorial_half*M/(double_fact*factorial);

        // calculate B
        for (int k = 0; k <=n; k++)
        {
            // calculate S0
            real factorial = 1;
            for (int i = 2; i <=n-k; i++)
                factorial *= i;

            real factorial2 = 1;
            for (int i = 2; i <= 2*n-k; i+=1)
                factorial2 *= i;

            real power2 = 1;
            for (int i = 0; i < n-k; i++)
                power2*=2;

            real S = factorial2/(factorial*power2);
            B[k] = S;

            // calculate B[k]
            for (int j = k; j < n; j++)
            {
                S *= (j+1)*(n-j)*2/((real)(j+1-k)*(2*n-j));
                B[k] += S;
            }
        }
    }

    InvertedKuzminToomreDisk::~InvertedKuzminToomreDisk()
    {
        delete[] B;
        delete[] P0;
        delete[] P1;
    }

    void InvertedKuzminToomreDisk::calculate_lambda_init(real const* y)
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

    void InvertedKuzminToomreDisk::calculate_lambda_run(real const* y)
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

    void InvertedKuzminToomreDisk::calculate_nu(const real* y)
    {
        real rho = y[RHO], z = y[Z];
        real abs_z = z>0? z:-z;
        real rb = sqrtl(rho*rho + (abs_z+b)*(abs_z+b));
        real P_arg = (abs_z + b)/rb;

        this->nu = 0;
        real b_pow = 1;
        real inv_rb = 1.0/rb;
        real inv_rb_pow = inv_rb;

        legendre_polynomials(P_arg, n+1, this->P0);

        for (int k = 0; k <= n; k++)
        {
            this->nu += B[k]*b_pow*inv_rb_pow*P0[k];
            b_pow *= -b;
            inv_rb_pow *= inv_rb;
        }
        this->nu *= N;
    }

    void InvertedKuzminToomreDisk::calculate_nu1(const real* y)
    {
        real rho = y[RHO], z = y[Z];
        real rho2 = rho*rho;
        int sign_z = z>0? 1:-1;
        real abs_z = sign_z*z;
        real rb = sqrtl(rho*rho + (abs_z+b)*(abs_z+b));
        real P_arg = (abs_z + b)/rb;

        nu = nu_rho = nu_z = 0;
        real b_pow = 1;
        real inv_rb = 1.0/rb;
        real inv_rb_pows = inv_rb;
        real inv_rb_pows2 = inv_rb*inv_rb*inv_rb;

        legendre_polynomials1(P_arg, n+1, this->P0, this->P1);

        for (int k = 0; k <= n; k++)
        {
            nu += B[k]*b_pow*inv_rb_pows*P0[k];
            nu_rho += B[k]*b_pow*inv_rb_pows2*rho*((k+1)*P0[k]+P_arg*P1[k]);
            nu_z += B[k]*b_pow*inv_rb_pows2*((k+1)*(abs_z+b)*P0[k]-rho2*inv_rb*P1[k]);

            b_pow *= -b;
            inv_rb_pows *= inv_rb;
            inv_rb_pows2 *= inv_rb;
        }
        nu *= N;
        nu_rho *= -N;
        nu_z *= -N*sign_z;
    }

    void InvertedKuzminToomreDisk::calculate_nu2(const real* y)
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

        // Calculate potential and derivatives
        real rho2 = rho*rho;
        int sign_z = z>0? 1:-1;
        real abs_z = sign_z*z;
        real rb = sqrtl(rho*rho + (abs_z+b)*(abs_z+b));
        real P_arg = (abs_z + b)/rb;

        nu = nu_rho = nu_z = 0;
        real b_pow = 1;
        real inv_rb = 1.0/rb;
        real inv_rb_pows = inv_rb;
        real inv_rb_pows2 = inv_rb*inv_rb*inv_rb;

        legendre_polynomials1(P_arg, n+1, this->P0, this->P1);

        for (int k = 0; k <= n; k++)
        {
            nu += B[k]*b_pow*inv_rb_pows*P0[k];
            nu_rho += B[k]*b_pow*inv_rb_pows2*rho*((k+1)*P0[k]+P_arg*P1[k]);
            nu_z += B[k]*b_pow*inv_rb_pows2*((k+1)*(abs_z+b)*P0[k]-rho2*inv_rb*P1[k]);

            b_pow *= -b;
            inv_rb_pows *= inv_rb;
            inv_rb_pows2 *= inv_rb;
        }
        nu *= N;
        nu_rho *= -N;
        nu_z *= -N*sign_z;
    }
}