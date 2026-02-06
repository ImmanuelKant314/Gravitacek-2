/**
 * @file stepperbase.hpp
 * @author Karel Kraus
 * @brief Class representing general stepper for solving ODEs.
 * 
 * @copyright Copyright (c) 2026
 */

#pragma once
#include "gravitacek2/setup.hpp"
#include "gravitacek2/integrator/odesystem.hpp"

#include <memory>

namespace gr2
{
    /**
     * @brief Base class for integrating ODEs.
     * 
     * Class enables us to calculate \f$\vec{y}(t + h)\f$ given \f$t\f$ and 
     * \f$\vec{y}(t)\f$. It is also possible to calculate error estimate of 
     * this value.
     */
    class StepperBase
    {
    protected:
        std::shared_ptr<OdeSystem> ode; //!<system of differential equations to be solved

        // ========== Values of y ==========
        real *y_in;     //!<input value of \f$\vec{y}\f$
        real *y_out;    //!<output value of \f$\vec{y}\f$
        real *y_err;    //!<estimates of error for each coordinate
        real *y_cur;    //!<current value of \f$\vec{y}\f$ during calculations
        real *y_help;   //!<value of \f$\vec{y}\f$ saved on the side

        // ========== Values of dydt ==========
        real *dydt_in;  //!<input value of \f$\vec{y}\f$
        real *dydt_out; //!<output value of \f$\vec{y}\f$
        real *dydt_cur; //!<array for storing current values of \f$\frac{\mathrm{d}\vec{y}}{\mathrm{d} t}\f$
        real *dydt_opt; //!<array for storing value of \f$\frac{\mathrm{d}\vec{y}}{\mathrm{d} t}\f$ for improving performance

        // ========== Internal constants ==========
        real t_in;  //!<value of time variable at the beginning of the step
        real h;     //!<length of the time step
        int n;      //!<number of differential equations
    public:
        /**
         * @brief Construct a new StepperBase object.
         * 
         */
        StepperBase();

        /**
         * @brief Destroy the StepperBase object.
         * 
         */
        virtual ~StepperBase();

        /**
         * @brief Set ODE for stepper.
         * 
         * Method saves ODE and allocates memory for the usage.
         * 
         * @param ode integrated OdeSystem
         */
        virtual void set_OdeSystem(std::shared_ptr<OdeSystem> ode);

        /**
         * @brief Calculate next step.
         * 
         * Method calculates values of \f$\vec{y}(t + h)\f$ given
         * \f$\vec{y}(t)\f$ and \f$t\f$. After calculation values are saved back
         * to `y`. 
         * 
         * If parameter `dydt_in` is given, value of derivative is used in the 
         * calculation. It is not neseccary, but it can reduce number of 
         * calculations.
         * 
         * If parametr `dydt_out` is given, value of derivative for time \f$t +
         * h\f$ saved there.
         * 
         * @param t time value
         * @param y coordinate values, also used to save final values
         * @param h time step
         * @param dense true if dense output should be calculated
         * @param dytd_in array for giving time derivative in the initial time (for accelerating calculations)
         * @param dydt_out array for storing time derivative in the final time (if given)
         */
        virtual void step(const real &t, real y[], const real &h, 
            const bool &dense=false, const real dydt_in[] = nullptr, 
            real dydt_out[] = nullptr) = 0;

        /**
         * @brief Calculate next step with error estimate.
         * 
         * Method calculates values of \f$\vec{y}(t + h)\f$ given
         * \f$\vec{y}(t)\f$ and \f$t\f$ with estimation of absolute error. After
         * calculation coordinate values are saved back to `y` and absolute
         * errors to `err`. If this method is not implemented, than error is
         * estimated by comparing one full step with two half steps. For
         * estimate of y two halfsteps are used.
         * 
         * If parameter `dydt_in` is given, value of derivative is used in the 
         * calculation. It is not neseccary, but it can reduce number of 
         * calculations.
         * 
         * If parametr `dydt_out` is given, value of derivative for time \f$t +
         * h\f$ saved there.
         * 
         * @note Default implementation is highly inefficient and should be used 
         * with caution. Also it can behave unexpectadly, so be careful with 
         * usage.
         * 
         * @param t time value
         * @param y coordinate values
         * @param h time step
         * @param err array to save estimates of absolute errors for each coordinate
         * @param dense true if dense output should be calculated
         * @param dytd_in array for giving time derivative in the initial time (for accelerating calculations)
         * @param dydt_out array for storing time derivative in the final time (if given)
         */
        virtual void step_err(const real &t, real y[], const real &h, 
            real err[], const bool& dense=false, const real dydt_in[] = nullptr, 
            real dydt_out[] = nullptr);

        /**
         * @brief Prepare stepper for dense output.
         */
        virtual void prepare_dense();

        /**
         * @brief Return dense output.
         * 
         * We use interpolation polynomial to calculate values.
         * 
         * @param i index of coordinate value
         * @param t time for evaluation
         * @return value of i-th coordinate given the time \f$t\f$
         */
        virtual real dense_out(const int& i, const real &t);

        /**
         * @brief Get order of precision for stepper.
         * 
         * Local error of stepper scales as
         * \f[
         * \text{error} \sim h^{\text{order} + 1}.
         * \f]
         * 
         * @return order of integration
         */
        virtual int get_order() const = 0;

        /**
         * @brief Get the order of error for stepper.
         * 
         * Local error of stepper scales as
         * \f[
         * \text{error} \sim h^{\text{error\_order}}.
         * \f]
         * 
         * @return order of error
         */
        virtual int get_err_order() const = 0;
    };
}