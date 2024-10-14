#pragma once
#include "gravitacek2/geomotion/geomotion.hpp"
#include "gravitacek2/geomotion/weyl.hpp"

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
        static const int UT = 4;        //!<index of four-velocity \f$u^t\f$
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

    /**
     * @brief Class representing geodesic motion in Schwarzschild spacetime in Weyl coordinates.
     * 
     * Spacetime is given by potential
     * \f[
     * \nu = \frac{1}{2} \log{\frac{d_1 + d_2 - 2M}{d_1 + d_2 + 2M}},
     * \f]
     * where \f$d_{1,2} = \sqrt{\rho^2 + \left(z \mp M\right)^2}\f$.
     * 
     * Metric function \f$\lambda\f$ is known exactly in this case and is equal to
     * \f[
     * \lambda = \frac{1}{2}\log{\frac{(d_1 + d_2)^2 - 4M^2}{d_1d_2}}.
     * \f]
     * 
     */
    class WeylSchwarzschild: public Weyl
    {
    protected:
        real M; //!<mass of the Swarzschild black hole

        virtual void calculate_lambda_exact(const real* y);
        virtual void calculate_lambda_integral(const real* y);
    public:
        /**
         * @brief Construct a new WeylSchwarzschild object.
         * 
         * @param M mass of black hole
         * @param lambda_exact true if value of \f$\lambda\f$ should be calculated exaxtly 
         */
        WeylSchwarzschild(real M = 1, LambdaEvaluation init=LambdaEvaluation::exact, LambdaEvaluation run=LambdaEvaluation::exact);

        virtual void calculate_lambda_init(const real* y) override;
        virtual void calculate_lambda_run(const real* y) override;

        virtual void calculate_nu(const real* y) override;
        virtual void calculate_nu1(const real* y) override;
        virtual void calculate_nu2(const real* y) override;
    };
}