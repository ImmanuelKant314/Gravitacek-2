#pragma once
#include "gravitacek2/setup.hpp"
#include "gravitacek2/odesolver/ode.hpp"

namespace gr2
{
    /**
     * @brief Base class for integrating ODE.
     * 
     * Class enables us to calculate \f$\vec{y}(t + h)\f$ given \f$t\f$ and 
     * \f$\vec{y}(t)\f$. It is also possible to calculate error estimate of 
     * this value.
     */
    class Stepper
    {
    protected:
        ODE *ode;
        real *yt, *yt2, *dydt, *dydt2;
        real t;
        int n;
    public:
        /**
         * @brief Construct a new Stepper object.
         * 
         */
        Stepper();

        /**
         * @brief Destroy the Stepper object.
         * 
         */
        ~Stepper();

        /**
         * @brief Set ODE for stepper.
         * 
         * Method saves ODE and allocates memory for the usage.
         * 
         * @param ode integrated ODE
         */
        virtual void set_ODE(ODE& ode);

        /**
         * @brief Reset stepper.
         * 
         * This function should be used after the next step is not continuation
         * of a previous one.
         * 
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
         * @param t time value
         * @param y coordinate values
         * @param h time step
         * @param err array to save estimates of absolute errors for each coordinate
         * @param dydt_in 
         * @param dydt_out
         */
        virtual void step(const real &t, real y[], const real &h, real err[], const real dydt_in[] = nullptr, real dydt_out[] = nullptr);

        /**
         * @brief Get the order of integration
         * 
         * If stepper was not used, it returns maximal order of integration. 
         * Else it return order of integration for the last step. The order can 
         * vary if the Stepper is adaptive.
         * 
         * @return int order of integration
         */
        virtual int get_order() const = 0;
    };
}