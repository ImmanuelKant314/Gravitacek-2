#pragma once
#include "gravitacek2/setup.hpp"
#include "gravitacek2/integrator/odesystem.hpp"

#include <memory>

namespace gr2
{
    /**
     * @brief Base class for integrating ODE.
     * 
     * Class enables us to calculate \f$\vec{y}(t + h)\f$ given \f$t\f$ and 
     * \f$\vec{y}(t)\f$. It is also possible to calculate error estimate of 
     * this value.
     */
    class StepperBase
    {
    protected:
        std::shared_ptr<OdeSystem> ode; //!<solved system of differential equations

        // ========== Values of y ==========
        real *y_in;     //!<input value of \f$\vec{y}\f$
        real *y_out;    //!<output value of \f$\vec{y}\f$
        real *y_err;    //!<estimate of \f$
        real *y_cur;    //!<current value of \f$\vec{y}\f$ during calculations
        real *y_help;   //!< value of \f$\vec{y}\f$ saved on the side

        // ========== Values of dydt ==========
        real *dydt_in;  //!<input value of \f$\vec{y}\f$
        real *dydt_out; //!<output value of \f$\vec{y}\f$
        real *dydt_cur; //!<array for storing current value of \f$\frac{\mathrm{d}\vec{y}}{\mathrm{d} t}\f$
        real *dydt_opt; //!<array for storing value of \f$\frac{\mathrm{d}\vec{y}}{\mathrm{d} t}\f$ for improving performance

        // ========== Internal constants ==========
        real t_in;      //!<value of time variable at the beginning of the step
        real h;         //!<lenght of the step
        int n;          //!<number of solved differential equations
        bool dense;     //!<calculation of dense output during steps

    public:
        /**
         * @brief Construct a new Stepper object.
         * 
         * @param dense Prepare dense output during integration.
         */
        StepperBase(bool dense);

        /**
         * @brief Destroy the Stepper object.
         * 
         */
        ~StepperBase();

        /**
         * @brief Set ODE for stepper.
         * 
         * Method saves ODE and allocates memory for the usage.
         * 
         * @param ode integrated ODE
         */
        virtual void set_OdeSystem(std::shared_ptr<OdeSystem> ode);

        /**
         * @brief Reset stepper.
         * 
         * This function should be used after the next step is not continuation
         * of a previous one.
         */
        virtual void reset()=0;

        /**
         * @brief Calculate next step.
         * 
         * Method calculates values of \f$\vec{y}\f(t + h)$ given $\f$\vec{y}(t)\f$ 
         * and \f$t\f$. After calculation values are saved back to `y`. 
         * 
         * If parameter `dydt_in` is given, value of derivative is used in the 
         * calculation. It is not neseccary, but it can reduce number of 
         * calculations.
         * 
         * If parametr `dydt_out` is given, value of derivative is returned.
         * 
         * @param t time value
         * @param y coordinate values, also used to save final values
         * @param h time step
         * @param dytd_in 
         * @param dydt_out
         */
        virtual void step(const real &t, real y[], const real &h, const real dydt_in[] = nullptr, real dydt_out[] = nullptr) = 0;

        /**
         * @brief Calculate next step with error estimate.
         * 
         * Method calculates values of \f$\vec{y}\f(t + h)$ given $\f$\vec{y}(t)\f$ 
         * and \f$t\f$ with estimation of absolute error. After calculation coordinate 
         * values are saved back to `y` and absolute errors to `err`. If this 
         * method is not implemented, than error is estimated by comparing one 
         * full step with two half steps. For estimate of y two halfsteps are 
         * used.
         * 
         * If parameter `dydt_in` is given, value of derivative is used in the 
         * calculation. It is not neseccary, but it can reduce number of 
         * calculations.
         * 
         * If parametr `dydt_out` is given, value of derivative is returned.
         * 
         * @note Default implementation is highly inefficient and should be used 
         * with caution. Also it can behave unexpectadly, so be careful with 
         * usage.
         * 
         * @param t time value
         * @param y coordinate values
         * @param h time step
         * @param err array to save estimates of absolute errors for each coordinate
         * @param dydt_in 
         * @param dydt_out
         */
        virtual void step_err(const real &t, real y[], const real &h, real err[], const real dydt_in[] = nullptr, real dydt_out[] = nullptr);

        /**
         * @brief Prepare stepper for dense output
         * 
         * @param h time step
         */
        virtual void prepare_dense();

        /**
         * @brief Dense output
         * 
         * @param i index of outputed value
         * @param h time step
         * @param theta fraction of time_step
         * @return vale of ith coordinate
         */
        virtual real dense_out(const int& i, const real &t);

        /**
         * @brief Get the order of integration
         * 
         * If stepper was not used, it returns maximal order of integration. 
         * Else it return order of integration for the last step. The order can 
         * vary if the Stepper is adaptive.
         * 
         * Local error of stepper scales as
         * \f[
         * \mbox{err} \sim h^{\mbox{order} + 1}.
         * \f]
         * 
         * @return order of integration
         */
        virtual int get_order() const = 0;

        /**
         * @brief Get the order of error.
         * 
         * Error is scaled as
         * \f[
         * \mbox{err} \sim h^\mbox{err_order}.
         * \f]
         * 
         * @return order of error
         */
        virtual int get_err_order() const = 0;
    };
}