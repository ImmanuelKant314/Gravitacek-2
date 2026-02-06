/**
 * @file event.hpp
 * @author Karel Kraus
 * @brief General event for integration.
 * 
 * @copyright Copyright (c) 2026
 * 
 */

#pragma once
#include "gravitacek2/setup.hpp"
#include "gravitacek2/integrator/stepperbase.hpp"

namespace gr2
{
    /**
     * @brief Type of event for integration.
     * 
     */
    enum EventType
    {
        data,       //!<record data
        modyfing,   //!<modify integration
    };

    /**
     * @brief Event for integrating ordinary differential equations.
     * 
     * When integrating ordinary differential equation it is necessary to record
     * data, find the intercection of trajectory with some set of point or to
     * stop integration when error is too high. For this purpose we use events.
     */
    class Event
    {
    protected:
        EventType type; //!<type of event
        bool terminal;  //!<is event terminal for integration
    public:
         /**
          * @brief Construct a new Event object.
          * 
          * @param type type of event
          * @param terminal true if event should stop the itnegration
          */
        Event(const EventType& type, const bool &terminal=false);

        /**
         * @brief Get type of event.
         * 
         */
        EventType get_type() const;
        
       /**
        * @brief Get if the event is terminal.
        * 
        */
        bool get_terminal() const;

        /**
         * @brief Return value of internal function of the event.
         * 
         * Each event uses some function for triggering. Event is triggered 
         * when the value of the internal function crosses 0.
         * 
         * @param t time variable
         * @param dt time step
         * @param y coordinate variables
         * @param dydt derivate of coordinate \f$\vec{y}\f$ with respect to \f$t\f$
         * @return value of internal function
         */
        virtual real value(const real &t, const real &dt, const real y[], const real dydt[]) = 0;

        /**
         * @brief Apply event.
         * 
         * When time and coordinates of event are found, event is applied. This
         * can save values or modify the integration of ordinary differential
         * equation.
         *  
         * @param stepper stepper to integrate ODE
         * @param t time variable
         * @param dt time step
         * @param y coordinate variables
         * @param dydt derivate of coordinate \f$\vec{y}\f$ with respect to \f$t\f$
         */
        virtual void apply(StepperBase* stepper, real &t, real &dt, real y[], real dydt[]) = 0;
    };
}