#pragma once
#include "gravitacek2/odesolver/stepper.hpp"

namespace gr2
{
    /**
     * @brief Stepper using algorithm RK4.
     * 
     */
    class RK4 : public Stepper
    {
    protected:
        real *k1;   // array of coefficients \f$k_1\f$ for Runge-Kutta method
        real *k2;   // array of coefficients \f$k_2\f$ for Runge-Kutta method
        real *k3;   // array of coefficients \f$k_3\f$ for Runge-Kutta method
        real *k4;   // array of coefficients \f$k_4\f$ for Runge-Kutta method
    public:
        RK4();
        ~RK4();
        virtual void set_ODE(ODE& ode);
        virtual void reset();
        virtual void step(const real &t, real y[], const real &h, const real dydt_in[] = nullptr, real dydt_out[] = nullptr) override;
        virtual int get_order() const override;
    };
}