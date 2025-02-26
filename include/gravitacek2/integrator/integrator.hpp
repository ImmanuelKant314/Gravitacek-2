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
    class Integrator
    {
    protected:
        // ========== Components of itnegrator ========== 
        OdeSystem *ode;                     //!<solved system of ordinary differential equations
        StepperBase *stepper;               //!<stepper used for integrating
        StepControllerBase *stepcontroller; //!<step controller

        // ========== Events ========== 
        std::vector<Event*> events_data;                //!<vector of data events
        std::vector<Event*> events_data_precise;        //!<vector of data events (precise)
        std::vector<Event*> events_modifying;           //!<vector of modyfying events
        std::vector<Event*> events_modifying_precise;   //!<vector of modyfying events
        std::vector<Event*> events_terminal;            //!<vector of terminal events
        std::vector<Event*> events_terminal_precise;    //!<vector of terminal events (precise)

        // ========== Values of events ========== 
        real* events_data_precise_values;               //!<tracked values of data events
        real* events_modifying_precise_values;          //!<tracked values of modifying events
        real* events_terminal_precise_values;           //!<tracked values of terminal events

        // ========== Number of events ========== 
        int number_of_events_data_precise;      //!<number of data events
        int number_of_events_modifying_precise; //!<number of modifying events (precise)
        int number_of_events_terminal_precise;  //!<number of terminal events (precise)

        // ========== Step size ========== 
        real h;         //!<current value of time step
        real h2;        //!<value of time step for trying new step
        real h3;        //!<value of time step for events

        // ========== Time variable ========== 
        real t;         //!<current value of time
        real t2;        //!<value of time for new time step
        real t3;        //!<value of time for events

        // ========== Coordinate variables ========== 
        real *yt;       //!<array for storing current value of \f$\vec{y}\f$
        real *yt2;      //!<array for trying next step
        real *yt3;      //!<array for trying events

        // ========== Derivative values ========== 
        real *dydt;     //!<array for calculating \f$\frac{\mathrm{d} \vec{y}}{\mathrm{d}}\f$
        real *dydt2;    //!<array for calculating \f$\frac{\mathrm{d} \vec{y}}{\mathrm{d}}\f$ when trying new step
        real *dydt3;    //!<array for calculating \f$\frac{\mathrm{d} \vec{y}}{\mathrm{d}}\f$ when calculating events

        // ========== Error values ========== 
        real *err;      //!<array for calculating error
        real *err2;     //!<array for calculating error when trying new step
        real *err3;     //!<array for calculating error in event

        /**
         * @brief Initializing basic variables.
         * 
         */
        void basic_setup();

        /**
         * @brief Ininialize stepper based on its name.
         * 
         * @param stepper_name Name of stepper
         */
        void init_stepper(const std::string& stepper_name);
        bool solve_event(Event *event, const real& previous_value_of_event);

    public:
        /**
         * @brief Construct a new Integrator object.
         * 
         * This constructor assumes constant steps in time.
         * 
         * @param ode set of ordinary differential equations to be solved
         * @param stepper_name name of stepper
         */
        Integrator(OdeSystem &ode, const std::string& stepper_name);

        /**
         * @brief Construct a new Integrator object.
         * 
         * This constructor assumes StandardStepController.
         * 
         * @param ode set of ordinary differential equations to be solved
         * @param stepper_name name of stepper
         * @param eps_abs expected absolute error
         * @param eps_rel expected relative error
         * @param a_y coefficient 
         * @param a_dydt coefficient
         */
        Integrator(std::shared_ptr<OdeSystem> ode, const std::string& stepper_name, const real &atol, const real &rtol);

        /**
         * @brief Destroy the Integrator object.
         * 
         */
        ~Integrator();

        /**
         * @brief Add event to simulation.
         * 
         * @param event 
         */
        void add_event(Event* event);

        /**
         * @brief Integrate ordinary differential equation
         * 
         * @param y_start 
         * @param t_start 
         * @param t_end 
         * @param h_start 
         */
        void integrate(const real y_start[], const real &t_start, const real &t_end, const real &h_start);
    };
}