#pragma once
#include "gravitacek2/geomotion/geomotion.hpp"

namespace gr2
{
    /**
     * @brief Class representing geodesic motion in axially symmetric 
     * Majumdar-Papapetrou spacetime. 
     * 
     * Spacetime is described in cylindrical coordinates by metric
     * \f[
     * g_{\mu\nu} = 
     * \begin{pmatrix}
     * -N^2 & 0 & 0 & 0 \\ 
     * 0 & \rho^2N^{-2} & 0 & 0 \\
     * 0 & 0 & N^{-2} & 0 \\
     * 0 & 0 & 0 & N^{-2}
     * \end{pmatrix}
     * \f]
     * \f$t\f$ is time coordinate and \f$\rho\f$, \f$\phi\f$, \f$z\f$ are coordinates
     * fo cylindrical type. \f$N\f$ is lapsa function and depends on \f$\rho\f$ and
     * \f$z\f$.
     */
    class MajumdarPapapetrouWeyl : public GeoMotion
    {
    protected:
        real N_inv;         //!<value of \f$N^{-1}\f$
        real N_inv_rho;     //!<value of \f$N^{-1}_{,\rho}\f$
        real N_inv_z;       //!<value of \f$N^{-1}_{,z}\f$
        real N_inv_rhorho;  //!<value of \f$N^{-1}_{,\rho\rho}\f$
        real N_inv_rhoz;    //!<value of \f$N^{-1}_{,\rho z}\f$
        real N_inv_zz;      //!<value of \f$N^{-1}_{,zz}\f$

    public:
        static const int T = 0;         //!<index of coordinate \f$t\f$
        static const int PHI = 1;       //!<index of coordinate \f$\phi\f$
        static const int RHO = 2;       //!<index of coordinate \f$\rho\f$
        static const int Z = 3;         //!<index of coordinate \f$z\f$
        static const int UT = 4;        //!<index of four-velocity \f$u^t\f$
        static const int UPHI = 5;      //!<index of four-velocity \f$u^\phi\f$
        static const int URHO = 6;      //!<index of four-velocity \f$u^\rho\f$
        static const int UZ = 7;        //!<index of four-velocity \f$u^z\f$

        /**
         * @brief Construct a new Majumdar Papapetrou Weyl object.
         * 
         */
        MajumdarPapapetrouWeyl();

        /**
         * @brief Destroy the Majumdar Papapetrou Weyl object.
         * 
         */
        virtual ~MajumdarPapapetrouWeyl();

        /**
         * @brief Calculate value of \f$N^{-1}\f$.
         * 
         * @param y coordinate values
         */
        virtual void calculate_N_inv(const real* y) = 0;

        /**
         * @brief Calculate value of \f$N^{-1}\f$ and its first derivatives.
         * 
         * @param y coordinate values
         */
        virtual void calculate_N_inv1(const real* y) = 0;

        /**
         * @brief Calculate value of \f$N^{-1}\f$ and its first and second derivatives.
         * 
         * @param y coordiante values
         */
        virtual void calculate_N_inv2(const real* y) = 0;

        /**
         * @brief Get value of \f$N^{-1}\f$.
         * 
         * @return real value of \f$N^{-1}\f$
         */
        real get_N_inv() const;

        /**
         * @brief Get value of \f$(N^{-1})_{,\rho}\f$.
         * 
         * @return real value of \f$(N^{-1})_{,\rho}\f$
         */
        real get_N_inv_rho() const;

        /**
         * @brief Get value of \f$(N^{-1}_{,z}\f$.
         * 
         * @return real value of \f$(N^{-1})_{,z}\f$
         */
        real get_N_inv_z() const;

        /**
         * @brief Get value of \f$(N^{-1})_{,\rho\rho}\f$.
         * 
         * @return real value of \f$(N^{-1})_{,\rho\rho}\f$
         */
        real get_N_inv_rhorho() const;

        /**
         * @brief Get value of \f$(N^{-1})_{,\rho z}\f$.
         * 
         * @return real value of \f$(N^{-1})_{,\rho z}\f$
         */
        real get_N_inv_rhoz() const;

        /**
         * @brief Get value of \f$(N^{-1})_{,zz}\f$.
         * 
         * @return real value of \f$(N^{-1})_{,zz}\f$
         */
        real get_N_inv_zz() const;

        // ========== Calculate tensors ========== 

        virtual void calculate_metric(const real *y) override;
        virtual void calculate_christoffel_symbols(const real *y) override;
        virtual void calculate_riemann_tensor(const real *y) override;
    };
}