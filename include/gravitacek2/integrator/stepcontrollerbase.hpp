#pragma once
#include "gravitacek2/setup.hpp"

namespace gr2
{
    /**
     * @brief Abstract class for changing step size.
     * 
     */
    class StepControllerBase
    {
    protected:
        int n; // number of solved ordinary differential equations
    public:
        /**
         * @brief Construct a new StepController object.
         * 
         * @param n number of equations
         */
        StepControllerBase(const int &n);

        /**
         * @brief Calculate new step size.
         * 
         * @param y coordinate values
         * @param err estimates for coordinate values
         * @param dydt derivatives of y with respect to t
         * @param h value of h, used as input and output
         * @return true if error is was small enought, false if it is necessary to repeat the step
         */
        virtual bool hadjust(const real y[], const real err[], const real dydt[], real &h) = 0;
    };
}