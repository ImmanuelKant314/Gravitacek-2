/**
 * @file integrator.hpp
 * @author Karel Kraus
 * @brief Integrator for solving ODEs.
 * 
 * @copyright Copyright (c) 2026
 */

#pragma once

#include "gravitacek2/setup.hpp"
#include "gravitacek2/integrator/odesystem.hpp"
#include "gravitacek2/integrator/stepperbase.hpp"
#include "gravitacek2/integrator/stepcontrollerbase.hpp"
#include "gravitacek2/integrator/event.hpp"

#include <vector>
#include <string>

namespace gr2
{
    /**
     * @brief Integrator for solving ODEs.
     * 
     * Integrator includes stepping algorithm, step-size controller and 
     * handling events.
     */
    class Integrator
    {
    protected:
        // ========== Components of integrator ========== 
        std::shared_ptr<OdeSystem> ode;     //!<solved system of ODEs
        StepperBase *stepper;               //!<stepper used for integration
        StepControllerBase *stepcontroller; //!<step controller

        // ========== Events ========== 
        std::vector<std::shared_ptr<Event>> events_data;        //!<vector of data events
        std::vector<std::shared_ptr<Event>> events_modifying;   //!<vector of modifying events
        
        // ========== Current event ==========
        std::shared_ptr<Event> current_event;   //!<event, which is being proceeded
        bool current_event_terminal;            //!<is the current event terminal?

        // ========== Values of events ========== 
        real* events_modifying_values;  //!<tracked values of modifying events

        // ========== Number of events ========== 
        int number_of_events_modifying; //!<number of modifying events (precise)

        // ========== Step size ========== 
        real h;     //!<current value of time step
        real h2;    //!<value of time step for trying new step
        real h3;    //!<value of time step for event

        // ========== Time variable ========== 
        real t;     //!<current value of time
        real t2;    //!<value of time for new time step
        real t3;    //!<value of time for event

        // ========== Coordinate variables ========== 
        real *yt;   //!<array for storing current value of \f$\vec{y}\f$
        real *yt2;  //!<array for trying next step
        real *yt3;  //!<array for trying events

        // ========== Derivative values ========== 
        real *dydt;     //!<array for calculating \f$\frac{\mathrm{d} \vec{y}}{\mathrm{d}}\f$
        real *dydt2;    //!<array for calculating \f$\frac{\mathrm{d} \vec{y}}{\mathrm{d}}\f$ when trying new step
        real *dydt3;    //!<array for calculating \f$\frac{\mathrm{d} \vec{y}}{\mathrm{d}}\f$ when calculating events

        // ========== Error values ========== 
        real *err;      //!<array for calculating error
        real *err2;     //!<array for calculating error when trying new step
        real *err3;     //!<array for calculating error in event

        // ========== Other variables ==========
        bool dense; //!<true if dense output should be used

        /**
         * @brief Initializing basic variables. 
         * 
         * This means creating new vectors and setting pointers to nullptr.
         */
        void basic_setup();

        /**
         * @brief Ininialize stepper based on its name.
         * 
         * @param stepper_name name of stepper
         */
        void init_stepper(const std::string& stepper_name);
        //TODO: This should be done using enum

        /**
         * @brief Try to trigger the event.
         * 
         * @param event event that we study
         * @param previous_value_of_event previous value of given event
         * @return true if event is triggered, else false
         */
        bool solve_event(std::shared_ptr<Event> event, const real& previous_value_of_event);

    public:
        /**
         * @brief Construct a new Integrator object.
         * 
         * This constructor assumes constant steps in time.
         * 
         * @param ode ODEs to be solved
         * @param stepper_name name of stepper
         */
        Integrator(std::shared_ptr<OdeSystem> ode, const std::string& stepper_name, const bool &dense=false);

        /**
         * @brief Construct a new Integrator object.
         * 
         * This constructor assumes StandardStepController.
         * 
         * @param ode ODEs to be solved
         * @param stepper_name name of stepper
         * @param eps_abs expected absolute error
         * @param eps_rel expected relative error
         * @param a_y coefficient for StandardStepController
         * @param a_dydt coefficient for StandardStepController
         */
        Integrator(std::shared_ptr<OdeSystem> ode, const std::string& stepper_name, const real &atol, const real &rtol, const bool &dense = false);

        /**
         * @brief Destroy the Integrator object.
         * 
         */
        ~Integrator();

        /**
         * @brief Add event for integration.
         * 
         * @param event event for integration
         */
        void add_event(std::shared_ptr<Event> event);

        /**
         * @brief Integrate ordinary differential equation
         * 
         * @param y_start initial coordinate values
         * @param t_start intial time
         * @param t_end final time
         * @param h_start initial time step
         */
        void integrate(const real y_start[], const real &t_start, const real &t_end, const real &h_start);
    };
}