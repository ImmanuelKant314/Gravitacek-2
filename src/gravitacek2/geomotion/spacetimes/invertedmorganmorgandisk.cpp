#include <stdexcept>
#include <cmath>

#include "gravitacek2/geomotion/spacetimes.hpp"
#include "gravitacek2/mymath.hpp"

namespace gr2
{
    void InvertedMorganMorganDisk::calculate_lambda_integral(const real* y)
    {
        calculate_lambda_from_inf_to_z(y, 1e-15);
    }

    InvertedMorganMorganDisk::InvertedMorganMorganDisk(int n, real M, real b, LambdaEvaluation init, LambdaEvaluation run): Weyl(init, run)
    {
        // save values
        this->n = n;
        this->M = M;
        this->b = b;

        // prepare arrays
        this->C = nullptr;
        this->P0 = nullptr;
        this->P1 = nullptr;
        this->Q0 = nullptr;
        this->Q1 = nullptr;

        // create arrays
        C = new real[n+1];
        P0 = new real[2*n+1];
        P1 = new real[2*n+1];
        Q0 = new real[n+1];
        Q1 = new real[n+1];

        // calculate normalization
        real pow2 = 2;
        real fact = 1;
        for (int i = 1; i<=n; i++)
        {
            pow2 *= 4;
            fact *= i;
        }
        N = -pow2*fact*fact*M/(pi*b);

        // calculate C[0]
        real fact2n1 = 1;
        for (int i = 2; i <= 2*n+1; i++)
            fact2n1 *= i;

        real CC = (real)1.0/(fact2n1);
        C[0] = CC;

        // calculate C[m]
        for (int m = 1; m <= n; m++)
        {
            CC *= -((real) 2*m*(2*m-1)*(n+m)*(n-m+1))/((real) (m*m*(2*n+2*m+1)*(2*n+2*m)));
            C[m] = (4*m+1)*CC;
        }
    }

    InvertedMorganMorganDisk::~InvertedMorganMorganDisk()
    {
        delete[] C, P0, P1, Q0, Q1;
    }

    void InvertedMorganMorganDisk::calculate_lambda_init(real const* y)
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

    void InvertedMorganMorganDisk::calculate_lambda_run(real const* y)
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

    void InvertedMorganMorganDisk::calculate_nu(const real *y)
    {
        real rho = y[RHO], z = y[Z];
        real alpha = rho*rho + z*z - b*b;
        real help_term = sqrtl(alpha*alpha + 4*(z*z)*(b*b));
        real x = sqrtl(alpha + help_term)/(sqrtl(2)*b);
        real y_ = sqrtl(-alpha + help_term)/(sqrtl(2)*b);
        help_term = sqrtl(x*x - y_*y_ + 1);
        real Y = y_/help_term;
        real X = x/help_term;

        legendre_polynomials(X, 2*n + 1, P0);
        special_function_Q2n(Y, n+1, Q0);

        nu = 0;
        for (int m = 0; m<= n; m++)
        {
            nu += C[m]*Q0[m]*P0[2*m];
        }
        nu *= N/help_term;
    }
    
    void InvertedMorganMorganDisk::calculate_nu1(const real *y)
    {
        real rho = y[RHO], z = y[Z];
        real alpha = rho*rho + z*z - b*b;
        real help_term1 = sqrtl(alpha*alpha + 4*(z*z)*(b*b));
        real x = sqrtl(alpha + help_term1)/(sqrtl(2)*b);
        real y_ = sqrtl(-alpha + help_term1)/(sqrtl(2)*b);
        real help_term2 = sqrtl(x*x - y_*y_ + 1);
        real Y = y_/help_term2;
        real X = x/help_term2;
        real x_rho = x*rho/help_term1;
        real y_rho = -y_*rho/help_term1;
        real x_z = z*(alpha + 2*b*b + help_term1)/(2*b*b*x*help_term1);
        real y_z = z*(alpha + 2*b*b - help_term1)/(2*b*b*y_*help_term1);
        real X_rho = ((1-y_*y_)*x_rho + x*y_*y_rho)/(help_term2*help_term2*help_term2);
        real Y_rho = (-x*y_*x_rho + (x*x+1)*y_rho)/(help_term2*help_term2*help_term2);
        real X_z = ((1-y_*y_)*x_z + x*y_*y_z)/(help_term2*help_term2*help_term2);
        real Y_z = (-x*y_*x_z + (x*x+1)*y_z)/(help_term2*help_term2*help_term2);

        legendre_polynomials1(X, 2*n + 1, P0, P1);
        special_function_Q2n1(Y, n + 1, Q0, Q1);

        nu = 0;
        nu_rho = 0;
        nu_z = 0;
        for (int m = 0; m<= n; m++)
        {
            nu += C[m]*Q0[m]*P0[2*m];

            nu_rho += C[m]*Q1[m]*P0[2*m]*Y_rho;
            nu_rho += C[m]*Q0[m]*P1[2*m]*X_rho;

            nu_z += C[m]*Q1[m]*P0[2*m]*Y_z;
            nu_z += C[m]*Q0[m]*P1[2*m]*X_z;
        }

        nu_rho += -rho*nu/(b*b)/(help_term2*help_term2);
        nu_z += -z*nu/(b*b)/(help_term2*help_term2);

        nu *= N/help_term2;
        nu_rho *= N/help_term2;
        nu_z *= N/help_term2;
    }

    void InvertedMorganMorganDisk::calculate_nu2(const real *y)
    {

    }
}

