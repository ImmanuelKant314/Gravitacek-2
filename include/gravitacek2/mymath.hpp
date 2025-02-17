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
    
    /**
     * @brief Calculate n Legendre polynomials.
     * 
     * @param x argument for Legendre polynomials
     * @param n number of calculated Legendre polynomials
     * @param p array for saving Legedre polynomials 
     */
    void legendre_polynomials(const real& x, const int& n, real* p);

    /**
    * @brief Calculate n Legendre polynomials with first derivatives.
    * 
    * @param x argument for Legendre polynomials
    * @param n number of calculated Legendre polynomials
    * @param p0 array for saving Legendre polynomials
    * @param p1 array for saving first derivatives of Legendre polynomials
    */
    void legendre_polynomials1(const real& x, const int& n, real* p0, real* p1);

    /**
     * @brief Calculate n values of special function \f$\mathcal{Q}_{2n}\f$.
     * 
     * Special function \f$\mathcal{Q}_n(x)\f$ is defined by the Legendre function of 
     * the second kind \f$Q_n(x)\f$ as
     * \f[
     * \mathcal{Q}_n(x) \equiv i Q_n(ix).
     * \f]
     * 
     * Important note, the definition of \f$Q_n\f$ is not typically used one, we use
     * \f[
     * Q_n(ix) = -i P_{2n}(ix) \left(\frac{\pi}{2} - \atan{x}\right)
     * - \sum_{j=1}^{2n}\frac{1}{j} P_{j-1}(ix) P_{2n-j}(ix).
     * \f]
     * 
     * @param x argument for special function
     * @param n number of calculated values
     * @param q array for saving values
     */
    void special_function_Q2n(const real& x, const int& n, real* q);

    /**
     * @brief Calculate n values of special function \f$\mathcal{Q}_{2n}\f$ and its derivatives.
     * 
     * Special function \f$\mathcal{Q}_n(x)\f$ is defined by the Legendre function of 
     * the second kind \f$Q_n(x)\f$ as
     * \f[
     * \mathcal{Q}_n(x) \equiv i Q_n(ix).
     * \f]
     * 
     * Important note, the definition of \f$Q_n\f$ is not typically used one, we use
     * \f[
     * Q_n(ix) = -i P_{2n}(ix) \left(\frac{\pi}{2} - \atan{x}\right)
     * - \sum_{j=1}^{2n}\frac{1}{j} P_{j-1}(ix) P_{2n-j}(ix).
     * \f]
     * 
     * @param x argument for special function
     * @param n number of calculated values
     * @param q0 array for saving values
     * @param q1 array for saving derivatives
     */
    void special_function_Q2n1(const real& x, const int& n, real* q0, real* q1);

    /**
     * @brief Numerically differentiate function using Richardson's extrapolation.
     * 
     * @tparam K order of extrapolation
     * @param func function to be differentiated
     * @param x point of differentiation
     * @param h0 initial step for differentiation
     * @param eps precision of differentiation
     * @return value of derivative
     */
    template<int K, class T>
    real richder(T func, const real& x, const real& h0, const real &eps = 10e-10)
    {
        const int JMIN = 5, JMAX = 20;
        int m;
        real R[K][JMAX];
        real h = 0.5*h0;
        R[0][0] = (func(x+h0) - func(x-h0))/(2*h0);
        for (int j = 1, n = 1; j < JMAX; j++, n *= 2, h*= 0.5)
        {
            // ========== Second order rule ==========
            R[0][j] = (func(x+h)-func(x-h))/(2*h);

            // ========== Richardsons's algorithm ==========
            m = 4;
            for (int i = 1; i < std::min(j, K); i++, m *= 4)
                R[i][j] = (m * R[i - 1][j] - R[i - 1][j - 1]) / (m - 1);

            // ========== Checking precision ==========
            if (j >= std::max(JMIN - 1, K + 1))
                if (fabsl(R[K - 1][j] - R[K - 1][j - 1]) < eps * (fabsl(R[K - 1][j - 1]) + 1))
                    return R[K - 1][j];

        }
        throw std::runtime_error("Too much iterations in routine richder");
    }

    /**
     * @brief Numerically calculate second derivative of the function using 
     * Richardson's extrapolation.
     * 
     * @tparam K order of extrapolation
     * @param func function to be differentiated
     * @param x point of differentiation
     * @param h0 initial step for differentiation
     * @param eps precision of differentiation
     * @return value of second derivative
     */
    template<int K, class T>
    real richder2(T func, const real& x, const real& h0, const real &eps = 10e-10)
    {
        const int JMIN = 5, JMAX = 20;
        int m;
        real R[K][JMAX];
        real h = 0.5*h0;
        real f0 = func(x);
        R[0][0] = (func(x+h0) -2*f0 + func(x-h0))/(h0*h0);
        for (int j = 1, n = 1; j < JMAX; j++, n *= 2, h*= 0.5)
        {
            // ========== Second order rule ==========
            R[0][j] = (func(x+h)-2*f0 + func(x-h))/(h*h);

            // ========== Richardsons's algorithm ==========
            m = 4;
            for (int i = 1; i < std::min(j, K); i++, m *= 4)
                R[i][j] = (m * R[i - 1][j] - R[i - 1][j - 1]) / (m - 1);

            // ========== Checking precision ==========
            if (j >= std::max(JMIN - 1, K + 1))
                if (fabsl(R[K - 1][j] - R[K - 1][j - 1]) < eps * (fabsl(R[K - 1][j - 1]) + 1))
                    return R[K - 1][j];

        }
        throw std::runtime_error("Too much iterations in routine richder2");
    }
}

