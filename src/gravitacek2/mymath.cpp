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
}