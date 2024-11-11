#include "gravitacek2/mymath.hpp"

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

}