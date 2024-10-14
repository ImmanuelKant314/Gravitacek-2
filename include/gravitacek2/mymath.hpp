#pragma once
#include "gravitacek2/setup.hpp"

#include <cmath>
#include <algorithm>
#include <exception>

namespace gr2
{
    /**
     * @brief Calculate complete elliptic integral of first and second kind.
     * 
     * Complete elliptic integral of the first kind is defined as
     * \f[
     * K(k) \equiv \int_0^{\pi/2} \frac{\dd\phi}{\sqrt{1 - k^2\sin^2{\phi}} }\dd \phi
     * \f]
     * 
     * Complete elliptic integral of the second kind is defined as 
     * \f[
     * E(k) \equiv \int_0^{\pi/2} \sqrt{1 - k^2\sin^2{\phi}}\dd \phi
     * \f]
     * 
     * @param k argument of elliptic integrals
     * @param K value of \f$K(k)\f$
     * @param E value of \f$E(k)\f$
     * @param eps precision
     */
    void elliptic_KE(const real& k, real &K, real &E, const real &eps = 1e-15);

    /**
     * @brief Integrate function using Romberg's integration.
     * 
     * @tparam K order of Romberg's method
     * @param func function to be integrated
     * @param a lower bound of integration
     * @param b upper bound of integration
     * @param eps presision of integration
     * @return value of integral
     */
    template<int K, class T>
    real romb(T func, const real &a, const real &b, const real &eps = 10e-10)
    {
        const int JMIN = 5, JMAX = 20;
        int m;
        real R[K][JMAX];
        real step = b - a, sum = 0, x;
        R[0][0] = 0.5 * step * (func(a) + func(b));
        for (int j = 1, n = 1; j < JMAX; j++, n *= 2, step *= 0.5)
        {
            // ========== Trapezoidal rule ==========
            sum = 0;
            x = a + 0.5 * step;
            for (int i = 0; i < n; i++, x += step)
                sum += func(x);
            R[0][j] = 0.5 * (R[0][j - 1] + sum * step);

            // ========== Romberg's algorithm ==========
            m = 4;
            for (int i = 1; i < std::min(j, K); i++, m *= 4)
                R[i][j] = (m * R[i - 1][j] - R[i - 1][j - 1]) / (m - 1);

            // ========== Checking precision ==========
            if (j >= std::max(JMIN - 1, K + 1))
                if (fabsl(R[K - 1][j] - R[K - 1][j - 1]) < eps * (fabsl(R[K - 1][j - 1]) + 1))
                    return R[K - 1][j];
        }
        throw std::runtime_error("Too much iterations in routine romb");
    }
}