/**
 * @file geomotion.hpp
 * @author Karel Kraus
 * @brief Class for representing ODE for geodesic motion in general space(-time).
 * 
 * @copyright Copyright (c) 2026
 */

#pragma once
#include "gravitacek2/integrator/odesystem.hpp"

namespace gr2
{
    /**
     * @brief General representation for ODEs describing geodesic motion.
     * 
     * Class has also additional functionality to calculate metric, Christoffel
     * symbols and Riemann tensor.
     */
    class GeoMotion : public OdeSystem
    {
    protected:
        int dim;                        //!<dimension of space(-time)

        //TODO: should be converted to 1D arrays
        real **metric;                  //!<metric tensor \f$g_{\mu\nu}\f$
        real ***christoffel_symbols;    //!<Christoffel symbols \f$\Gamma^{\mu}_{\kappa\lambda}\f$
        real ****riemann_tensor;        //!<Riemann tensor \f$R^{\mu}_{\nu\kappa\lambda}\f$

        real *y_m;                      //!<position where `metric` is calculated
        real *y_c;                      //!<position where `christoffel_symbols`is calculated
        real *y_r;                      //!<position where `riemann_tensor`is calculated

        /**
         * @brief Check if calculation in given space(-time) point is necessary.
         * 
         * @param y current coordinate variables
         * @param y_save coordinate variables, where caclulation was done previously
         * @param n dimension of `y` and `y_save`
         * @return true if calculation is necessary
         * @return false if calculation is not necessary
         */
        bool necessary_calculate(const real *y, real *&y_save, const int& n);
    public:
        // ========== Constructors & destructors ========== 
        /**
         * @brief Construct a new GeoMotion object.
         * 
         * @param dim dimension of space(-time)
         * @param n number of differential equations, typically 2*`dim`
         */
        GeoMotion(const int &dim, const int &n);

        /**
         * @brief Destroy the GeoMotion object.
         * 
         */
        virtual ~GeoMotion();

        // ========== Calculate ========== 

        /**
         * @brief Calculate metric.
         * 
         * Metric is calculated is the form
         * \f$
         * g_{\mu\nu}.
         * \f$
         * 
         * @param y coordinate variables
         */
        virtual void calculate_metric(const real *y) = 0;

        /**
         * @brief Calculate Christoffel symbols.
         * 
         * Christoffel symbols are calculated in the form
         * \f[
         * \tensor{\Gamma}{^\mu_\kappa_\lambda} = 
         * \frac{1}{2}g^{\mu\sigma}\left(
         * g_{\sigma\kappa,\lambda} + g_{\lambda\sigma,\kappa} - g_{\kappa\lambda,\sigma}
         * \right).
         * \f]
         * 
         * @param y coordinate variables
         */
        virtual void calculate_christoffel_symbols(const real *y) = 0;

        /**
         * @brief Calculate Riemann tensor.
         * 
         * Riemann tensor is calculated in the form
         * \f[
         * \tensor{R}{^\mu_\nu_\kappa_\lambda} = 
         * \tensor{\Gamma}{^\mu_\nu_\lambda_{,\kappa}} -
         * \tensor{\Gamma}{^\mu_\nu_\kappa_{,\lambda}} +
         * \tensor{\Gamma}{^\mu_\rho_\kappa} \tensor{\Gamma}{^\rho_\nu_\lambda} -
         * \tensor{\Gamma}{^\mu_\rho_\lambda} \tensor{\Gamma}{^\rho_\nu_\kappa}.
         * \f]
         * 
         * @param y coordinate variables
         */
        virtual void calculate_riemann_tensor(const real *y) = 0;
        
        // ========== Getters ========== 

        /**
         * @brief Get spatial dimensions of space(-time).
         * 
         * @return dimension
         */
        int get_dim() const;

        /**
         * @brief Get metric tensor.
         * 
         * Metric tensor is in form \f$g_{\mu\nu}\f$.
         * 
         * @return metric tensor
         */
        real **get_metric() const;

        /**
         * @brief Get Christoffel symbols.
         * 
         * Christoffel symbols are in form \f$\tensor{\Gamma}{^\mu_\nu_\kappa}\f$.
         * 
         * @return Christoffel symbols
         */
        real ***get_christoffel_symbols() const;

        /**
         * @brief Get Riemann tensor.
         * 
         * Riemann tensor is in form of \f$\tensor{R}{^\mu_\nu_\kappa_\lambda}\f$.
         * 
         * @return Riemann tensor
         */
        real ****get_riemann_tensor() const;

        // ========== Function ========== 

        /**
         * @brief Calculate time derivative of state vector.
         * 
         * State vector \f$\vec{y}\f$ contain coordinate values and velocities:
         * \f[
         * \vec{y} = (x^\mu, u^\mu).
         * \f]
         * Derivative of state vector is calculated as
         * \f[
         * \dv{\vec{y}}{t} = \left(u^\mu, -\tensor{\Gamma}{^\mu_\kappa_\lambda}u^\kappa u^\lambda\right).
         * \f]
         * 
         * @param t time variable
         * @param y state vector
         * @param dydt derivation of state vector with respect to \f$t\f$
         */
        virtual void function(const real &t, const real y[], real dydt[]) override;
    };
}