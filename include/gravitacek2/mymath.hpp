#include "gravitacek2/setup.hpp"

#include <cmath>

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
}