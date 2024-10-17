// ========== include - my library ========== 
#include "gravitacek2/integrator/integrator.hpp"
#include "gravitacek2/integrator/steppers.hpp"
#include "gravitacek2/integrator/stepcontrollers.hpp"

// ========== include - standard libraries ========== 
#include <stdexcept>
#include <cmath>

// ========== macros ========== 
#define MAX_ITERATIONS_HADJUST 10
#define MAX_ITERATIONS_SOLVE_EVENT 20
#define EVENT_PRECISION 1e-15

namespace gr2
{
    void Integrator::basic_setup()
    {
        this->ode = nullptr;
        this->stepper = nullptr;
        this->stepcontroller = nullptr;

        this->events_data_precise_values = nullptr;
        this->events_modifying_precise_values = nullptr;
        this->events_terminal_precise_values = nullptr;

        this->yt = nullptr;
        this->yt2 = nullptr;
        this->yt3 = nullptr;

        this->dydt = nullptr;
        this->dydt2 = nullptr;
        this->dydt3 = nullptr;

        this->err = nullptr;
        this->err2 = nullptr;
        this->err3 = nullptr;

        this->events_data = std::vector<Event*>();
        this->events_modifying = std::vector<Event*>();
        this->events_modifying_precise = std::vector<Event*>();
        this->events_terminal = std::vector<Event*>();
        this->events_terminal_precise = std::vector<Event*>();
    }

    void Integrator::init_stepper(const std::string& stepper_name)
    {
        delete stepper;
        if (stepper_name == "RK4")
            this->stepper = new RK4();
        else
            throw std::invalid_argument("no integrator with given name found");
    }

    bool Integrator::solve_event(Event *event, const real& previous_value_of_event)
    {
        // prepare values
        int i;

        // coppy array yt2 to yt3 and h2 to h3
        for (i = 0; i< this->ode->get_n(); i++)
        {
            yt3[i] = yt2[i];
        }
        h3 = h2;
        t3 = t2;

        // calculate current value of event
        real current_value_of_event = event->value(t3, yt3, dydt3);

        // check if event is triggered
        // if event is not triggered return true
        // else if current value of event is close enough, return true
        if (current_value_of_event*previous_value_of_event > 0 || current_value_of_event != 0)
            return false;
        else if (fabs(previous_value_of_event) < EVENT_PRECISION) 
            return true;

        // prepare values for secant method
        real a = previous_value_of_event, b = current_value_of_event;
        real h_a = 0, h_b = h3;

        // ========== Secant method ==========
        for (i = 0; i < MAX_ITERATIONS_SOLVE_EVENT; i++)
        {
            h3 = (h_a*b-h_b*a)/(b-a); // new value of step size
            t3 = t + h3;
            // copy starting value of yt to yt3
            for (int j = 0; j < this->ode->get_n(); j++)
                yt3[j] = yt[j];

            // take step and calculate new value of event
            this->stepper->step_err(t, yt3, h3, err3, dydt, dydt3);
            current_value_of_event = event->value(t3, yt3, dydt3);

            // reduce interval for h
            if (current_value_of_event*previous_value_of_event > 0)
            {
                h_a = h3;
                a = current_value_of_event;
            }
            else
            {
                h_b = h3;
                b = current_value_of_event;
                // if new value of event is close enought, stop for-cycle
                if( (h_b-h_a) < EVENT_PRECISION*std::max(h_a, h_b))
                    return true;
            }
        }
       
        if (i == MAX_ITERATIONS_SOLVE_EVENT)
            throw std::runtime_error("precise time of event could not be found");

        return true;
    }
        
    Integrator::Integrator(OdeSystem &ode, const std::string& stepper_name)
    {
        this->basic_setup();
        this->ode = &ode;
        this->init_stepper(stepper_name);
        this->stepper->set_OdeSystem(ode);
    }

    Integrator::Integrator(OdeSystem &ode, const std::string& stepper_name, const real &atol, const real &rtol)
    {
        this->basic_setup();
        this->ode = &ode;
        this->init_stepper(stepper_name);
        this->stepper->set_OdeSystem(ode);
        delete this->stepcontroller;
        this->stepcontroller = new StepControllerNR(ode.get_n(), this->stepper->get_err_order(), atol, rtol);
    }

    void Integrator::add_event(Event* event)
    {
        switch (event->get_type())
        {
        case EventType::data:
            this->events_data.push_back(event);
            break;
        case EventType::data_precise:
            this->events_data_precise.push_back(event);
            break;
        case EventType::modyfing:
            this->events_modifying.push_back(event);
            break;
        case EventType::modyfing_precise:
            this->events_modifying_precise.push_back(event);
            break;
        case EventType::terminal:
            this->events_terminal.push_back(event);
            break;
        case EventType::terminal_precise:
            this->events_terminal_precise.push_back(event);
            break;
        default:
            throw std::invalid_argument("invalid type of event");
            break;
        }
    }

    Integrator::~Integrator()
    {
        delete stepper;
        delete stepcontroller;
    }

    void Integrator::integrate(const real y_start[], const real &t_start, const real &t_end, const real &h_start)
    {
        // prepare variables
        int i;
        Event *current_event = nullptr;
        int n = this->ode->get_n();

        this->yt = new real[n];
        this->yt2 = new real[n];
        this->yt3 = new real[n];
        this->dydt = new real[n];
        this->dydt2 = new real[n];
        this->dydt3 = new real[n];
        this->err = new real[n];
        this->err2 = new real[n];
        this->err3 = new real[n];

        // copy values internaly
        for (int i = 0; i < n; i++)
            this->yt[i] = this->yt2[i] = y_start[i];
        this->t = t_start;
        this->h = this->h2 = this->h3 = h_start;
        t2 = t3 = t + h2;

        // number of events
        number_of_events_data_precise = events_data.size();
        number_of_events_modifying_precise = events_modifying.size();
        number_of_events_terminal_precise = events_terminal.size();

        // prepare event values
        events_data_precise_values = new real[number_of_events_data_precise];
        events_modifying_precise_values = new real[number_of_events_modifying_precise];
        events_terminal_precise_values = new real[number_of_events_terminal_precise];

        // prepare values of events
        for (int i = 0; i < number_of_events_data_precise; i++)
            events_data_precise_values[i] = events_data[i]->value(t_start, yt, dydt);
        for (int i = 0; i < number_of_events_modifying_precise; i++)
            events_modifying_precise_values[i] = events_modifying[i]->value(t_start, yt, dydt);
        for (int i = 0; i < number_of_events_terminal_precise; i++)
            events_terminal_precise_values[i] = events_terminal[i]->value(t_start, yt, dydt);

        // cycle for calculating new values of y
        while (t < t_end)
        {   
            // taky a step
            this->stepper->step_err(t, yt2, h, err2);

            // null current event
            current_event = nullptr;

            // check modifying events
            for (int i = 0; i < number_of_events_modifying_precise; i++)
            {
                if(this->solve_event(events_modifying[i], events_modifying_precise_values[i]))
                {
                    for (int i = 0; i < n; i++)
                    {
                        yt2[i] = yt3[i];
                    }
                    h2 = h3;
                    t2 = t3;
                    current_event = events_modifying[i];
                }
            }

            // check terminal events
            for (int i = 0; i < number_of_events_terminal_precise; i++)
            {
                if(this->solve_event(events_terminal[i], events_terminal_precise_values[i]))
                {
                    for (int i = 0; i < this->ode->get_n(); i++)
                    {
                        yt2[i] = yt3[i];
                    }
                    h2 = h3;
                    t2 = t3;
                    current_event = events_terminal[i];
                }
            }

            i = 0;
            if (this->stepcontroller)
            {
                for (i = 0; i < MAX_ITERATIONS_HADJUST; i++)
                {
                    if(this->stepcontroller->hadjust(this->yt2, this->err2, this->dydt, this->h2))
                        break;
                    for (int j = 0; j < n; j++)
                        yt2[j] = yt[j];
                    this->stepper->step_err(t, yt2, h2, err2, dydt, dydt2);
                    t2 = t + h2;
                }
            }

            // if too much iteration for hadjust, throw exception
            if (i >= MAX_ITERATIONS_HADJUST)
                throw std::runtime_error("optimal step size was not found, MAX_ITERATIONS_HADJUST reached");

            for (int i = 0; i < number_of_events_data_precise; i++)
            {
                if(this->solve_event(events_data[i], events_data_precise_values[i]))
                    events_data[i]->apply(t3, yt3, dydt3);
            }

            // "commit" to the step
            for (int i = 0; i < n; i++)
            {
                yt[i] = yt2[i];
                dydt[i] = dydt2[i];
            }
            t = t2;
            if (this->stepcontroller)
                h = h2;
            h2 = h;
            h3 = h;
            t2 = t3 = t + h;

            // apply event
            if (current_event)
                current_event->apply(t, yt, dydt);
            
            // not precise events
            for (auto &event : events_data)
                if (event->value(t, yt, dydt) == 0)
                    event->apply(t, yt, dydt);

            for (auto &event : events_modifying)
                if (event->value(t, yt, dydt) == 0)
                    event->apply(t, yt, dydt);

            bool terminal = false;
            for(auto &event : events_terminal)
                if (event->value(t, yt, dydt)==0)
                {
                    event->apply(t, yt, dydt);
                    terminal = true;
                }

            // kill precise event
            if (current_event && current_event->get_type() == EventType::terminal)
                break;

            // kill unprecise event
            if (terminal)
                break;

            // calculate new values of events
            for (int i = 0; i < number_of_events_data_precise; i++)
                events_data_precise_values[i] = events_data[i]->value(t, yt, dydt);
            for (int i = 0; i < number_of_events_modifying_precise; i++)
                events_modifying_precise_values[i] = events_modifying[i]->value(t, yt, dydt);
            for (int i = 0; i < number_of_events_terminal_precise; i++)
                events_terminal_precise_values[i] = events_terminal[i]->value(t, yt, dydt);

            
        }

        // delete events
        delete[] events_data_precise_values, events_modifying_precise_values, events_terminal_precise_values;
        // delete other arrays
        delete[] yt, yt2, yt3, dydt, dydt2, dydt3, err, err2, err3;
    }
}