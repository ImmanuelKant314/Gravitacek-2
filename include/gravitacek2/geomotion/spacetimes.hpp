#pragma once
#include "gravitacek2/geomotion/geomotion.hpp"
#include "gravitacek2/geomotion/weyl.hpp"
#include "gravitacek2/geomotion/majumadpapapetrouweyl.hpp"

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
        WeylSchwarzschild(real M = 1, LambdaEvaluation init=LambdaEvaluation::exact, LambdaEvaluation run=LambdaEvaluation::diff);

        virtual void calculate_lambda_init(const real* y) override;
        virtual void calculate_lambda_run(const real* y) override;

        virtual void calculate_nu(const real* y) override;
        virtual void calculate_nu1(const real* y) override;
        virtual void calculate_nu2(const real* y) override;
    };

    /**
     * @brief Class representing geodesic motion in spacetime of Bach-Weyl ring.
     * 
     * Spacetime is given by potential
     * \f[
     * \nu = -\frac{2\mathcal{M} K(k)}{\pi l_2},
     * \f]
     * where \f$K(k) \equiv \int_0^{\pi/2}\frac{\dd \phi}{\sqrt{1-k^2\sin^2{\phi}}}\f$
     * is complete elliptic integral of the first kind, \f$l_{1,2} = \sqrt{(\rho \mp b)^2 + z^2)}\f$
     * and \f$k = \frac{\sqrt{4\rho b}}{l_2}\f$. \f$\mathcal{M}\f$ is mass of the ring 
     * and \f$b\f$ is a radius of the ring.
     * 
     */
    class BachWeylRing: public Weyl
    {
    protected:
        real M; //!<mass of the Bach-Weyl ring
        real b; //!<radius of the Bach-Weyl ring

        virtual void calculate_lambda_integral(const real* y);
    public:
        BachWeylRing(real M, real b, LambdaEvaluation init=LambdaEvaluation::integral, LambdaEvaluation run=LambdaEvaluation::diff);

        virtual void calculate_lambda_init(const real* y) override;
        virtual void calculate_lambda_run(const real* y) override;

        virtual void calculate_nu(const real* y) override;
        virtual void calculate_nu1(const real* y) override;
        virtual void calculate_nu2(const real* y) override;
    };

    /**
     * @brief Class representing geodesic motion in space-time of inverted Kuzmin-Toomre disks.
     * 
     * Spacetime is given by potential
     * \f[
     * \nu_\text{iKT}^{(n)} = -\binom{n+1/2}{n} \frac{\mathcal{M}}{(1+2n)!!}\sum_{k=0}^{n} \mathcal{B}^{(n)}_k \frac{(-b)^k}{r_b^{k+1}}P_k\left(\frac{|z|+b}{r_b}\right),
     * \f]
     * where
     * \f[
     * r_b^2 \equiv \rho^2 + (|z| + b)^2
     * \f]
     * and
     * \f[
     * \mathcal{B}_k^{(n)} \equiv \sum_{j=k}^n \binom{j}{k}\frac{(2n-j)!}{2^{n-j}(n-j)!}.
     * \f]
     * \f$n\f$ is index of the disk (in the family), \f$\mathcal{M}\f$ 
     * represents mass of the disk, \f$b\f$ refers to dypical radial scale of the disk. 
     * \f$P_k\f$ refers to Legendre polynomial.
     */
    class InvertedKuzminToomreDisk: public Weyl
    {
    protected:
        int n;  //!<index of inverted Kuzmin-Toomre disk
        real M; //!<mass of the inverted Kuzmin-Toomre disk
        real b; //!<radius of the inverted Kuzmin-Toomre disk 

        real N;     //!<normalization constant
        real *B;    //!<coefficients for potential
        real *P0;   //!<values of Legendre polynomials
        real *P1;   //!<values of derivatives of Legendre polynomials

        virtual void calculate_lambda_integral(const real* y);
    public:
        InvertedKuzminToomreDisk(int n, real M, real b, LambdaEvaluation init=LambdaEvaluation::integral, LambdaEvaluation run=LambdaEvaluation::diff);
        ~InvertedKuzminToomreDisk();

        virtual void calculate_lambda_init(const real* y) override;
        virtual void calculate_lambda_run(const real* y) override;

        virtual void calculate_nu(const real* y) override;
        virtual void calculate_nu1(const real* y) override;
        virtual void calculate_nu2(const real* y) override;
    };

    class InvertedMorganMorganDisk: public Weyl
    {
    protected:
        int n;  //!<index of inverted Kuzmin-Toomre disk
        real M; //!<mass of the inverted Kuzmin-Toomre disk
        real b; //!<radius of the inverted Kuzmin-Toomre disk 

        real N;     //!<normalization constant
        real *C;    //!<constants for calculating potential
        real *P0;   //!<values of Legendre polynomials
        real *P1;   //!<values of derivatives of Legendre polynomials
        real *Q0;   //!<values of special function Q
        real *Q1;   //!<values of derivatives of special function Q

        virtual void calculate_lambda_integral(const real* y);
    public:
        InvertedMorganMorganDisk(int n, real M, real b, LambdaEvaluation init=LambdaEvaluation::integral, LambdaEvaluation run=LambdaEvaluation::diff);
        ~InvertedMorganMorganDisk();

        virtual void calculate_lambda_init(const real* y) override;
        virtual void calculate_lambda_run(const real* y) override;

        virtual void calculate_nu(const real* y) override;
        virtual void calculate_nu1(const real* y) override;
        virtual void calculate_nu2(const real* y) override;
    };

    class ReissnerNordstromMPW : public MajumdarPapapetrouWeyl
    {
    protected:
        real M;
    public:
        ReissnerNordstromMPW(const real& M);
        ~ReissnerNordstromMPW();

        virtual void calculate_N_inv(const real* y) override;
        virtual void calculate_N_inv1(const real* y) override;
        virtual void calculate_N_inv2(const real* y) override;
    };

    class MajumdarPapapetrouRing : public MajumdarPapapetrouWeyl
    {
    protected:
        real M;
        real b;
    public:
        MajumdarPapapetrouRing(const real& M, const real &b);
        ~MajumdarPapapetrouRing();

        virtual void calculate_N_inv(const real* y) override;
        virtual void calculate_N_inv1(const real* y) override;
        virtual void calculate_N_inv2(const real* y) override;
    };
}