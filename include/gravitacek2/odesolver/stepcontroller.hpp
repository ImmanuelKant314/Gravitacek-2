#pragma once
#include "gravitacek2/setup.hpp"

namespace gr2
{
    /**
     * @brief Abstract class for changing stepsize.
     * 
     */
    class StepController
    {
    protected:
        int n; // number of ODEs
    public:
        /**
         * @brief Construct a new StepController object.
         * 
         * @param n number of equations
         */
        StepController(const int &n);

        /**
         * @brief Calculate new step size.
         * 
         * @param y coordinate values
         * @param err estimates for coordinate values
         * @param dydt derivatives of y with respect to t
         * @param h previous step size
         * @return real new step size
         */
        virtual real hadjust(const real y[], const real err[], const real dydt[], const real &h) = 0;
    };
}