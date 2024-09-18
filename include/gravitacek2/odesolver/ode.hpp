#pragma once
#include "gravitacek2/setup.hpp"

namespace gr2
{
    /**
     * @brief Class representing Ordinary Differential Equations (ODE).
     *
     * Assume we have equation 
     * \f[
     * \frac{\mathrm{d}\vec{y}}{\mathrm{d}t} = f(t, \vec{y}),
     * \f]
     * where \f$t\f$ is a time variable and \f$\vec{y}\f$ are coordinate 
     * variables. This class represents the right side of the equation. 
     * It enables us to calculate derivation of \f$\vec{y}\f$ wth respect to 
     * \f$t\f$ with given \f$t\f$ and \f$\vec{y}\f$.
     * 
     * The class is abstract and is used to create concrete differential 
     * equations.
     */
    class ODE
    {
    protected:
        int n;
    public:
        /**
         * @brief Construct a new ODE object
         * 
         * @param n number of equations
         */
        ODE(const int &n);

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
         * @param dydt array for returning values of derivation of \f$\vec{y}\f$ with respect to \f$t\f$
         */
        virtual void function(const gr2::REAL &t, const gr2::REAL y[], gr2::REAL dydt[]) = 0;
    };
}