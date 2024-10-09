#pragma once
#include "gravitacek2/geomotion/geomotion.hpp"

namespace gr2
{
    /**
     * @brief Class representing geodesic motion in Schwarzschild spacetime.
     * 
     * Spacetime is represented by metric
     * \f[
     * g_{\mu\nu} = 
     * \begin{pmatrix}
     * -1+\frac{2M}{r} & 0 & 0 & 0 \\ 
     * 0 & \frac{1}{1-\frac{2M}{r}} & 0 & 0 \\ 
     * 0 & 0 & r^2 & 0 \\
     * 0 & 0 & 0 & r^2 \sin^2\theta
     * \end{pmatrix},
     * \f]
     * where \f$t\f$ is time coordinate and \f$r\f$, \f$\theta\f$, \f$\phi\f$ 
     * are coordinates of spherical type.
     */
    class Schwarzschild : public GeoMotion
    {
    protected:
        real M;                         //!<mass of the Swarzschild black hole
    public:
        static const int T = 0;         //!<index of coordinate \f$t\f$
        static const int R = 1;         //!<index of coordinate \f$r\f$
        static const int THETA = 2;     //!<index of coordinate \f$\theta\f$
        static const int PHI = 3;       //!<index of coordinate \f$\varphi\f$
        static const int UT = 4;        //!<index of four-velocit \f$u^t\f$
        static const int UR = 5;        //!<index of four-velocity \f$u^r\f$
        static const int UTHETA = 6;    //!<index of four-velocity \f$u^\theta\f$
        static const int UPHI = 7;      //!<index of four-velocity \f$u^\varphi\f$

        /**
         * @brief Construct a new Schwarzschild object
         * 
         * @param M mass of black hole
         */
        Schwarzschild(real M = 1);

        virtual void calculate_metric(const real *y) override;
        virtual void calculate_christoffel_symbols(const real *y) override;
        virtual void calculate_riemann_tensor(const real *y) override;
    };
}