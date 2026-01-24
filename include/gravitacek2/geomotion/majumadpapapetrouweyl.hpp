/**
 * @file majumadpapapetrouweyl.hpp
 * @author Karel Kraus
 * @brief Class representing general Majumdar-Papapetrou spacetime in Weyl
 * coordinates.
 * 
 * @copyright Copyright (c) 2026
 */

#pragma once
#include "gravitacek2/geomotion/geomotion.hpp"
#include <vector>
#include <memory>

namespace gr2
{
    /**
     * @brief GeoMotion class for general axially symmetric Majumdar-Papapetrou
     * space-time in Weyl coordiantes.
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
     * \f$t\f$ is time coordinate and \f$\rho\f$, \f$\phi\f$, \f$z\f$ are
     * coordinates of cylindrical type. \f$N\f$ is lapse function and depends on
     * \f$\rho\f$ and \f$z\f$.
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
         * @brief Calculate value of \f$N^{-1}\f$ and its first and second
         * derivatives.
         * 
         * @param y coordiante values
         */
        virtual void calculate_N_inv2(const real* y) = 0;

        /**
         * @brief Get value of \f$N^{-1}\f$.
         * 
         * @return value of \f$N^{-1}\f$
         */
        real get_N_inv() const;

        /**
         * @brief Get value of \f$(N^{-1})_{,\rho}\f$.
         * 
         * @return value of \f$(N^{-1})_{,\rho}\f$
         */
        real get_N_inv_rho() const;

        /**
         * @brief Get value of \f$(N^{-1})_{,z}\f$.
         * 
         * @return value of \f$(N^{-1})_{,z}\f$
         */
        real get_N_inv_z() const;

        /**
         * @brief Get value of \f$(N^{-1})_{,\rho\rho}\f$.
         * 
         * @return value of \f$(N^{-1})_{,\rho\rho}\f$
         */
        real get_N_inv_rhorho() const;

        /**
         * @brief Get value of \f$(N^{-1})_{,\rho z}\f$.
         * 
         * @return value of \f$(N^{-1})_{,\rho z}\f$
         */
        real get_N_inv_rhoz() const;

        /**
         * @brief Get value of \f$(N^{-1})_{,zz}\f$.
         * 
         * @return value of \f$(N^{-1})_{,zz}\f$
         */
        real get_N_inv_zz() const;

        // ========== Calculate tensors ========== 

        virtual void calculate_metric(const real *y) override;
        virtual void calculate_christoffel_symbols(const real *y) override;
        virtual void calculate_riemann_tensor(const real *y) override;
    };

    /**
     * @brief GeoMotion class for superposition of axially symmetric
     * Majumdar-Papapetrou space-times in Weyl coordinates.
     * 
     * For lapse function \f$N\f$ of \f$n\f$ individual Majumdar-Papapetrou
     * space-times it holds that
     * \f[
     * \frac{1}{N} = \sum_{i=1}^n \frac{1}{N_i} - n + 1,
     * \f]
     * where $N_i$ are lapse functions for individual space-times.
     */
    class CombinedMPW : public MajumdarPapapetrouWeyl
    {
    protected:
        std::vector<std::shared_ptr<MajumdarPapapetrouWeyl>> sources;   //!<vector of individual spacetimes
    public:
        /**
         * @brief Construct a new CombinedMPW object.
         * 
         * @param sources vector of individual sources
         */
        CombinedMPW(std::vector<std::shared_ptr<MajumdarPapapetrouWeyl>> sources);
        ~CombinedMPW();

        virtual void calculate_N_inv(const real* y) override;
        virtual void calculate_N_inv1(const real* y) override;
        virtual void calculate_N_inv2(const real* y) override;
    };
}