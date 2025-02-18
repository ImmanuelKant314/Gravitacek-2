#pragma once
#include "gravitacek2/geomotion/geomotion.hpp"

namespace gr2
{
    enum LambdaEvaluation
    {
        exact,
        diff,
        integral,
        custom1,
        custom2,
        custom3
    };

    /**
     * @brief Class representing geodesic motion in general Weyl spacetime.
     * 
     * Spacetime is described by metric
     * \f[
     * g_{\mu\nu} = 
     * \begin{pmatrix}
     * - e^{2\nu} & 0 & 0 & 0 \\
     * 0 & \rho^2e^{-2\nu} & 0 & 0 \\
     * 0 & 0 & e^{2\lambda - 2\nu} & 0 \\
     * 0 & 0 & 0 & e^{2\lambda - 2\nu}
     * \end{pmatrix}.
     * \f]
     * \f$t\f$ is time coordinate and \f$\rho\f$, \f$\phi\f$, \f$z\f$ are coordinates
     * of cylindrical type. \f$\nu\f$ and \f$\lambda\f$ are  metric functions 
     * dependent on \f$\rho\f$ and \f$z\f$.
     * 
     * Metric function \f$\nu\f$ satisfies Laplace equation in form of
     * \f[
     * \Delta \nu = \frac{1}{\rho}\nu_{,\rho} + \nu_{,\rho\rho} + \nu_{,zz} = 0.
     * \f]
     * Metric function \f$\lambda\f$ can be described using derivatives
     * \f[
     * \begin{align*}
     * \lambda_{,\rho} &= \rho \left[ (\nu_{,\rho})^2 - (\nu_{,z})^2\right],\\
     * \lambda_{,z} &= 2\rho \nu_{,\rho}\nu_{,z}.
     * \end{align*}
     * \f]
     */
    class Weyl : public GeoMotion
    {
    protected:
        LambdaEvaluation lambda_eval_init;   //!<type of calculation of \f$\lambda\f$ when initializing
        LambdaEvaluation lambda_eval_run;    //!<type of calculation of \f$\lambda\f$ when running simulation
        
        int lambda_index;               //!<index where value of \f$\lambda\f$ should be stored

        real nu;            //!<value of \f$\nu\f$
        real nu_rho;        //!<value of \f$\nu_{,\rho}\f$
        real nu_z;          //!<value of \f$\nu_{,z}\f$
        real nu_rhorho;     //!<value of \f$\nu_{,\rho\rho}\f$
        real nu_rhoz;       //!<value of \f$\nu{,\rho z}=\nu{,z\rho}\f$
        real nu_zz;         //!<value of \f$\nu{,zz}\f$

        real lambda;        //!<value of \f$\lambda\f$
        real lambda_rho;    //!<value of \f$\lambda_{,\rho}\f$
        real lambda_z;      //!<value of \f$\lambda_{,z}\f$
        real lambda_rhorho; //!<value of \f$\lambda_{,\rho\rho}\f$
        real lambda_rhoz;   //!<value of \f$\lambda_{,\rho z}\f$
        real lambda_zz;     //!<value of \f$\lambda_{,zz}\f$

        /**
         * @brief Calculate value of \f$\lambda\f$ by integrating from \f$z = \inf\f$ to \f$z = z_0\f$.
         * 
         * Integration is done by Romberg's algorithm with K=5. Function integrated
         * is
         * \f[
         * \lambda(\rho, z) = -\int_{0}^{\frac{1}{1+|z|}} \frac{1}{t^2}\lambda_{,z}\left(\rho, \frac{1}{t} - 1\right) \dd t.
         * \f]
         * 
         * @param y coordinate values
         */
        void calculate_lambda_from_inf_to_z(const real* y, const real& eps=1e-10);
        void calculate_lambda_diff(const real* y);

    public:
        static const int T = 0;         //!<index of coordinate \f$t\f$
        static const int PHI = 1;       //!<index of coordinate \f$\phi\f$
        static const int RHO = 2;       //!<index of coordinate \f$\rho\f$
        static const int Z = 3;         //!<index of coordinate \f$z\f$
        static const int UT = 4;        //!<index of four-velocity \f$u^t\f$
        static const int UPHI = 5;      //!<index of four-velocity \f$u^\phi\f$
        static const int URHO = 6;      //!<index of four-velocity \f$u^\rho\f$
        static const int UZ = 7;        //!<index of four-velocity \f$u^z\f$
        static const int LAMBDA = 8;    //!<index of value \f$\lambda\f$

        /**
         * @brief Construct a new Weyl object.
         * 
         * @param lambda_exact true if value of \f$\lambda\f$ should be calculated exaxtly
         */
        Weyl(LambdaEvaluation init, LambdaEvaluation run);

        /**
         * @brief Destroy the Weyl object.
         * 
         */
        virtual ~Weyl();

        /**
         * @brief Calculate value of \f$\nu\f$.
         * 
         * @param y coordinate values
         */
        virtual void calculate_nu(const real* y) = 0;

        /**
         * @brief Calculate value of \f$\nu\f$ and its first derivatives.
         * 
         * @param y coordinate values
         */
        virtual void calculate_nu1(const real *y) = 0;

        /**
         * @brief Calculate value of \f$\nu\f$ and its first and second derivatives.
         * 
         * @param y coordinate values
         */
        virtual void calculate_nu2(const real *y) = 0;

        /**
         * @brief Get value of \f$\nu\f$.
         * 
         * @return real value of \f$\nu\f$
         */
        real get_nu() const;

        /**
         * @brief Get value of \f$\nu_{,\rho}\f$.
         * 
         * @return real value of \f$\nu_{,\rho}\f$
         */
        real get_nu_rho() const;

        /**
         * @brief Get value of \f$\nu_{,z}\f$.
         * 
         * @return real value of \f$\nu_{,z}\f$
         */
        real get_nu_z() const;

        /**
         * @brief Get value of \f$\nu_{,\rho\rho}\f$.
         * 
         * @return real value of \f$\nu_{,\rho\rho}\f$
         */
        real get_nu_rhorho() const;

        /**
         * @brief Get value of \f$\nu_{,\rho z}\f$.
         * 
         * @return real value of \f$\nu_{,\rho z}\f$
         */
        real get_nu_rhoz() const;

        /**
         * @brief Get value of \f$\nu_{,zz}\f$.
         * 
         * @return real value of \f$\nu_{,zz}\f$
         */
        real get_nu_zz() const;

        /**
         * @brief Calculate value of \f$\lambda\f$.
         * 
         * Value of \f$\lambda\f$ is usually calculated by integral or exactly if
         * formula is known.
         * 
         * @param y coordinate variables
         */
        virtual void calculate_lambda_init(const real* y) = 0;

        /**
         * @brief Calculate value of \f$\lambda\f$ optimally.
         * 
         * If `lambda_exact = true` then value is calculated exactly, else
         * value is extracted from `y`.
         * 
         * @param y  coordinate variables (may include value of \f$\lambda\f$)
         */
        virtual void calculate_lambda_run(const real* y) = 0;

        /**
         * @brief Get value of \f$\lambda\f$
         * 
         * @return real value of \f$\lambda\f$
         */
        real get_lambda() const;

        /**
         * @brief Set index, where value of \f$\lambda\f$ is saved.
         * 
         * When calculating with `y` and `lambda_exact = false`, value of 
         * \f$\lambda\f$ has to be given using `y[lambda_index]`.
         * 
         * @param lambda_index index of value \f$\lambda\f$
         */
        void set_lambda_index(const int& lambda_index);

        /**
         * @brief Get index, where value of \f$\lambda\f$ is saved.
         * 
         * When calculating with `y` and `lambda_exact = false`, value of 
         * \f$\lambda\f$ has to be given using `y[lambda_index]`.
         * 
         * @return value of lambda_index 
         */
        int get_lambda_index() const;

        // ========== Calculate tensors ========== 

        virtual void calculate_metric(const real *y) override;
        virtual void calculate_christoffel_symbols(const real *y) override;
        virtual void calculate_riemann_tensor(const real *y) override;

    };
}