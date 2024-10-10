#pragma once
#include "gravitacek2/odesolver/ode.hpp"

namespace gr2
{
    /**
     * @brief Class GeoMotion represents ODE for geodesic motion.
     * 
     */
    class GeoMotion : public ODE
    {
    protected:
        int dim;                        //!<dimension of space

        real **metric;                  //!<metric tensor \f$g_{\mu\nu}\f$
        real ***christoffel_symbols;    //!<Christoffel symbols \f$\Gamma^{\mu}_{\kappa\lambda}\f$
        real ****riemann_tensor;        //!<Riemann tensor \f$R^{\mu}_{\nu\kappa\lambda}\f$

        real *y_m;                      //!<position where `metric` is calculated
        real *y_c;                      //!<position where `christoffel_symbols`is calculated
        real *y_r;                      //!<position where `riemann_tensor`is calculated

        /**
         * @brief Check if calculation in this point is necessary.
         * 
         * @param y current coordinate variables
         * @param y_save coordinate variables, where caclulation was done previously
         * @param n dimension of `y` and `y_save`
         * @return true if calculation is necessary
         * @return false if calculation is not necessary
         */
        bool necesarry_calculate(const real *y, real *y_save, const int& n);
    public:
        // ========== Constructors & destructors ========== 
        /**
         * @brief Construct a new GeoMotion object.
         * 
         * @param dim dimension of space
         * @param n number of differential equations, typically 2*`dim`
         */
        GeoMotion(const int &dim, const int &n);

        /**
         * @brief Destroy the GeoMotion object.
         * 
         */
        ~GeoMotion();

        // ========== Calculate ========== 

        /**
         * @brief Calculate metric.
         * 
         * Metric is calculated is form
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
         * Christoffel symbols are calculated in form
         * \f$
         * \Gamma^{\mu}_{\nu\lambda}.
         * \f$
         * 
         * @param y coordinate variables
         */
        virtual void calculate_christoffel_symbols(const real *y) = 0;

        /**
         * @brief Calculate Riemann tensor.
         * 
         * Riemann tensor is calculated in form
         * \f$
         * R^{\mu}_{\nu\kappa\lambda}.
         * \f$
         * 
         * @param y coordinate variables
         */
        virtual void calculate_riemann_tensor(const real *y) = 0;
        
        // ========== Getters ========== 

        /**
         * @brief Get metric tensor.
         * 
         * Metric tensor is in form \f$g_{\mu\nu}\f$.
         * 
         * @return real** metric tensor
         */
        real **get_metric() const;

        /**
         * @brief Get Christoffel symbols.
         * 
         * Christoffel symbols are in form \f$\Gamma^{\mu}_{\nu\kappa}\f$.
         * 
         * @return real*** Christoffel symbols
         */
        real ***get_christoffel_symbols() const;

        /**
         * @brief Get Riemann tensor.
         * 
         * Riemann tensor is in form of \f$R^{\mu}_{\nu\kappa\lambda}\f$.
         * 
         * @return real**** Riemann tensor
         */
        real ****get_riemann_tensor() const;

        // ========== Function ========== 

        virtual void function(const real &t, const real y[], real dydt[]) override;
    };
}