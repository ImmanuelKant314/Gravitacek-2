#pragma once
#include "gravitacek2/setup.hpp"

namespace gr2
{
    /**
     * @brief Represent type of event.
     * 
     */
    enum EventType
    {
        data,
        data_precise,
        modyfing,
        modyfing_precise,
        terminal,
        terminal_precise
    };

    /**
     * @brief Event for integrating ordinary differential equations.
     * 
     * When integrating ordinary differential equation it is necessary to record 
     * data, find the intercection of trajectory with some set of point or to stop
     * integration when error is too high. For this reason there is this class.
     */
    class Event
    {
    protected:
        EventType type; //!<type of event
    public:
        /**
         * @brief Construct a new Event object.
         * 
         * @param type type of event
         */
        Event(const EventType& type);

        /**
         * @brief Get type of Event object.
         * 
         */
        EventType get_type() const;

        /**
         * @brief Return value of event.
         * 
         * Each event has a continuous function. Event is triggerd if function 
         * crosses 0.
         * 
         * @param t time variable
         * @param y coordinate variables
         * @param err error of current step
         * @param dydt derivate of coordinate \f$\vec{y}\f$ with respect to \f$t\f$
         * @param h step size
         * @return value of internal function
         */
        virtual real value(const real &t, const real y[], const real dydt[]) = 0;

        /**
         * @brief Apply event.
         * 
         * When time and coordinates of event are found, event is applied. This 
         * can save values or alter integration of ordinary differential 
         * equation.
         *  
         * @param t time variable
         * @param y coordinate variables
         * @param err error of current step
         * @param dydt derivate of coordinate \f$\vec{y}\f$ with respect to \f$t\f$
         * @param h step size
         * @return value of internal function
         */
        virtual void apply(real &t, real y[], real dydt[]) = 0;
    };
}