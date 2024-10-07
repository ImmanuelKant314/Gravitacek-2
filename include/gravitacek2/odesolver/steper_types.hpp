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
        virtual int get_err_order() const override;
    };

    /**
     * @brief Stepper using algorithm DoPr853
     * 
     */
    class DoPr853 : public Stepper
    {
    protected:
        real *k1, *k2, *k3, *k4, *k5, *k6, *k7, *k8, *k9, *k10, *k11, *k12;
        real *k_help;
    public:
        DoPr853();
        ~DoPr853();
        virtual void set_ODE(ODE& ode);
        virtual void reset();
        virtual void step(const real &t, real y[], const real &h, const real dydt_in[] = nullptr, real dydt_out[] = nullptr) override;
        virtual void step_err(const real &t, real y[], const real &h, real err[], const real dydt_in[] = nullptr, real dydt_out[] = nullptr) override;
        virtual int get_order() const override;
        virtual int get_err_order() const override;
    };
}