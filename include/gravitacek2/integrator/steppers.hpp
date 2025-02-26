#pragma once
#include "gravitacek2/integrator/stepperbase.hpp"

namespace gr2
{
    /**
     * @brief Stepper using algorithm RK4.
     * 
     */
    class RK4 : public StepperBase
    {
    protected:
        real *k1;   //!<array of coefficients \f$k_1\f$ for Runge-Kutta method
        real *k2;   //!<array of coefficients \f$k_2\f$ for Runge-Kutta method
        real *k3;   //!<array of coefficients \f$k_3\f$ for Runge-Kutta method
        real *k4;   //!<array of coefficients \f$k_4\f$ for Runge-Kutta method
    public:
        /**
         * @brief Construct a new RK4 object.
         * 
         */
        RK4();

        /**
         * @brief Destroy the RK4 object.
         * 
         */
        ~RK4();

        virtual void set_OdeSystem(std::shared_ptr<OdeSystem> ode) override;
        virtual void reset();
        virtual void step(const real &t, real y[], const real &h, const bool& dense=false, const real dydt_in[] = nullptr, real dydt_out[] = nullptr) override;
        virtual int get_order() const override;
        virtual int get_err_order() const override;
    };

    /**
     * @brief Stepper using algorithm DoPr853
     * 
     */
    class DoPr853 : public StepperBase
    {
    protected:
        real *k1;   //!<array of coefficients \f$k_1\f$ for Runge-Kutta method
        real *k2;   //!<array of coefficients \f$k_2\f$ for Runge-Kutta method
        real *k3;   //!<array of coefficients \f$k_3\f$ for Runge-Kutta method
        real *k4;   //!<array of coefficients \f$k_4\f$ for Runge-Kutta method 
        real *k5;   //!<array of coefficients \f$k_5\f$ for Runge-Kutta method
        real *k6;   //!<array of coefficients \f$k_6\f$ for Runge-Kutta method 
        real *k7;   //!<array of coefficients \f$k_7\f$ for Runge-Kutta method
        real *k8;   //!<array of coefficients \f$k_8\f$ for Runge-Kutta method
        real *k9;   //!<array of coefficients \f$k_9\f$ for Runge-Kutta method
        real *k10;  //!<array of coefficients \f$k_{10}\f$ for Runge-Kutta method 
        real *k11;  //!<array of coefficients \f$k_{11}\f$ for Runge-Kutta method
        real *k12;  //!<array of coefficients \f$k_{12}\f$ for Runge-Kutta method
        real *k_help;
    public:
        /**
         * @brief Construct a new DoPr853 object.
         * 
         */
        DoPr853();

        /**
         * @brief Destroy the DoPr853 object.
         * 
         */
        ~DoPr853();

        virtual void set_OdeSystem(std::shared_ptr<OdeSystem> ode) override;
        virtual void reset();
        virtual void step(const real &t, real y[], const real &h, const bool &dense=false, const real dydt_in[] = nullptr, real dydt_out[] = nullptr) override;
        virtual void step_err(const real &t, real y[], const real &h, real err[], const bool &dense=false, const real dydt_in[] = nullptr, real dydt_out[] = nullptr) override;
        virtual int get_order() const override;
        virtual int get_err_order() const override;
    };
}