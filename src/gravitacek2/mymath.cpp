#include "gravitacek2/mymath.hpp"
#include <cmath>

namespace gr2
{
    void elliptic_KE(const real& k, real &K, real &E, const real &eps)
    { 
        real x = sqrtl(1 - k * k), y = 1;
        real x0 = x, y0 = y;
        real newx , newy; 
        real fac = 0.5, sum = 0; 
        real eeps = 2.7 * sqrtl(eps); 
        while (fabsl(x - y) > eeps * fabsl(x)) 
        {  
            newx = 0.5 * (x + y); 
            newy = sqrtl(x * y); 
            x = newx; 
            y = newy; 
            sum += fac * (x - y) * (x - y); 
            fac *= 2; 
        }  
        K = pi / (x + y);
        E = ((x0 + y0) * (x0 + y0) / 4 - sum) * K; 
    }

    void legendre_polynomials(const real& x, const int& n, real*p)
    {
        // evalueate \f$P_0(x)\f$
        if (n > 0)
            p[0] = 1;
        
        // evalueate \f$P_1(x)\f$
        if (n > 1)
            p[1] = x;

        // apply recurent formula
        for (int i = 1; i < n-1; i++)
            p[i+1] = ((2*i+1)*x*p[i] - i*p[i-1])/(i+1);
    }

    void legendre_polynomials1(const real& x, const int& n, real *p0, real *p1)
    {
        // evaluate \f$P_0(x)\f$ and \f$P_0'(x)\f$
        if (n > 0)
        {
            p0[0] = 1;
            p1[0] = 0;
        }
        
        // evaluate \f$P_1(x)\f$ and \f$P_1'(x)\f$
        if (n > 1)
        {
            p0[1] = x;
            p1[1] = 1;
        }

        // apply recurent formulas
        for (int i = 1; i < n-1; i++)
        {
            p0[i+1] = ((2*i+1)*x*p0[i] - i*p0[i-1])/(i+1);
            p1[i+1] = (i+1)*p0[i] + x*p1[i];
        }
    }

    void special_function_Q2n(const real& x, const int& n, real* q)
    {
        gr2::real q0, q1;
        if (n > 1)
        {
            q0 = pi_2-std::atan(x);
            q[0] = q0;
        }
        q1 = x*q0-1;
        for (int i = 1, j = 2; i<n; i++, j++)
        {
            q0 = (-(2*j-1)*x*q1-(j-1)*q0)/((real) j);
            q[i] = q0;
            j++;
            q1 = ((2*j-1)*x*q0-(j-1)*q1)/j;
        }
    }

    void special_function_Q2n1(const real& x, const int& n, real* q0, real* q1)
    {
        gr2::real q0_val, q1_val;
        gr2::real q0_der, q1_der;

        if (n > 1)
        {
            q0_val = pi_2-std::atan(x);
            q0_der = -1.0/(x*x+1);
            q0[0] = q0_val;
            q1[0] = q0_der;
        }

        if (n > 2)
        {
            q1_val = x*q0_val-1;
            q1_der = q0_val + x*q0_der;
        }

        for (int i = 1, j = 2; i<n; i++, j++)
        {
            q0_val = (-(2*j-1)*x*q1_val-(j-1)*q0_val)/((real) j);
            q0[i] = q0_val;
            q0_der = (-(2*j-1)*(x*q1_der + q1_val) - (j-1)*q0_der)/((real) j);
            q1[i] = q0_der;
            j++;
            q1_val = ((2*j-1)*x*q0_val-(j-1)*q1_val)/j;
            q1_der = ((2*j-1)*(x*q0_der + q0_val)-(j-1)*q1_der)/j;
        }
    }

}