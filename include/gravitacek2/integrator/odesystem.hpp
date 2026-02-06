/**
 * @file odesystem.hpp
 * @author Karel Kraus
 * @brief Classes for representing Ordinary Differential Equations.
 * 
 * @copyright Copyright (c) 2026
 */

#pragma once
#include "gravitacek2/setup.hpp"
#include <vector>
#include <memory>

namespace gr2
{
    /**
     * @brief Represenation of general Ordinary Differential Equations (ODEs).
     *
     * Assume we have equation 
     * \f[
     * \frac{\mathrm{d}\vec{y}}{\mathrm{d}t} = f(t, \vec{y}),
     * \f]
     * where \f$t\f$ is a time variable and \f$\vec{y}\f$ are coordinate
     * variables. This class represents the right side of the equation. It
     * enables us to calculate derivative of \f$\vec{y}\f$ with respect to
     * \f$t\f$ with values \f$t\f$ and \f$\vec{y}\f$.
     */
    class OdeSystem
    {
    protected:
        int n; //!<number of ordinary differential equations
    public:
        /**
         * @brief Construct a new OdeSystem object
         * 
         * @param n number of equations
         */
        OdeSystem(const int &n);

        virtual ~OdeSystem();

        /**
         * @brief Return number of equations.
         *
         * @return number of equations
         */
        int get_n() const;

        /**
         * @brief Calculate derivative of \f$\vec{y}\f$ with respect to \f$t\f$.
         * 
         * After calculation value is saved to `dydt`.
         * 
         * @param t time value
         * @param y coordinate values
         * @param dydt array for returning values of derivative of \f$\vec{y}\f$
         * with respect to \f$t\f$
         */
        virtual void function(const real &t, const real y[], real dydt[]) = 0;
    };

    /**
     * @brief OdeSystem constructed as a combination of more OdeSystems.
     * 
     * This class can combine mode OdeSystems which can be used in two different ways.
     * Firstly we can solve more separate ODEs at once, for example we can
     * simulate trajectory of two different particles. Secondly, we can results
     * from one ODE to calculate the second one, for example we can calculate 
     * trajectory for one particle and for the second one we can use
     * linearization (assuming the second particle will be always close to the 
     * first one).
     * 
     */
    class CombinedOdeSystem : public OdeSystem
    {
    protected:
        std::vector<std::shared_ptr<OdeSystem>> odes; //!< OdeSystems to be solved together
    public:

        /**
         * @brief Construct a new CombinedOdeSystem object.
         * 
         * @param odes ODEs to be solved together
         */
        CombinedOdeSystem(std::vector<std::shared_ptr<OdeSystem>> odes);
        void function(const real &t, const real y[], real dydt[]) override;
    };
}