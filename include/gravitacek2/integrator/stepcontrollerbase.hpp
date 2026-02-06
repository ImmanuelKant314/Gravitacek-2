/**
 * @file stepcontrollerbase.hpp
 * @author Karel Kraus
 * @brief Class for managing stepsize during integration of ODEs.
 * 
 * @copyright Copyright (c) 2026
 */

#pragma once
#include "gravitacek2/setup.hpp"

namespace gr2
{
    /**
     * @brief Abstract class for changing step size during integration of ODEs.
     * 
     */
    class StepControllerBase
    {
    protected:
        int n; //!<number of solved ordinary differential equations
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
         * @param err estimates of error for coordinate values
         * @param dydt derivatives of \f$\vec{y}\f$ with respect to \f$t\f$
         * @param h value of stepsize, used as input and output
         * @return true if error is small enought, false if it is necessary to repeat the step
         */
        virtual bool hadjust(const real y[], const real err[], const real dydt[], real &h) = 0;
    };
}