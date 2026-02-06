/**
 * @file steppers.hpp
 * @author Karel Kraus
 * @brief Specific steppers for solving ODEs.
 * 
 * @copyright Copyright (c) 2026
 */

#pragma once
#include "gravitacek2/integrator/stepperbase.hpp"

namespace gr2
{
    /**
     * @brief Stepper using algorithm RK4.
     * 
     * Coefficient for this stepper can be found in the source code or for 
     * example in the wikipedia article <a
     * href=https://en.wikipedia.org/wiki/Runge%E2%80%93Kutta_methods#Examples>Runge–Kutta
     * methods</a>.
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
     * @brief Stepper using algorithm DoPr853.
     * 
     * Coefficient for the algorithm can be found in the source code or in the 
     * <a href="https://www.unige.ch/~hairer/software.html">original
     * implementation</a>.
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

        real *pc1;  //!<array of coefficients `pc1` for dense output
        real *pc2;  //!<array of coefficients `pc2` for dense output
        real *pc3;  //!<array of coefficients `pc3` for dense output
        real *pc4;  //!<array of coefficients `pc4` for dense output
        real *pc5;  //!<array of coefficients `pc5` for dense output
        real *pc6;  //!<array of coefficients `pc6` for dense output
        real *pc7;  //!<array of coefficients `pc7` for dense output
        real *pc8;  //!<array of coefficients `pc8` for dense output

        real *k_help;   //!<array of coefficients for calculating error
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
        virtual void prepare_dense() override;
        virtual real dense_out(const int& i, const real &t) override;
    };
}