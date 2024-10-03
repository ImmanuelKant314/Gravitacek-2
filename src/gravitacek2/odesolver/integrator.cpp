// ========== include - my library ========== 
#include "gravitacek2/odesolver/integrator.hpp"
#include "gravitacek2/odesolver/steper_types.hpp"
#include "gravitacek2/odesolver/stepcontroller_types.hpp"

// ========== include - standard libraries ========== 
#include <stdexcept>
#include <cmath>

// ========== macros ========== 
#define MAX_ITERATIONS_HADJUST 20
#define MAX_ITERATIONS_SOLVE_EVENT 20
#define EVENT_PRECISION 1e-12
#define BIAS 1.0

namespace gr2
{
    void Integrator::basic_setup()
    {
        this->ode = nullptr;
        this->stepper = nullptr;
        this->stepcontroller = nullptr;

        this->events_data_values = nullptr;
        this->events_modifying_values = nullptr;
        this->events_terminal_values = nullptr;

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
        this->events_terminal = std::vector<Event*>();
    }

    void Integrator::init_stepper(const std::string& stepper_name)
    {
        delete stepper;
        if (stepper_name == "RK4")
            this->stepper = new RK4();
        else
            throw std::invalid_argument("no integrator with given name found");
    }

    bool Integrator::solve_event(Event *event, real& previous_value_of_event)
    {
        int i;
        for (i = 0; i< this->ode->get_n(); i++)
        {
            yt3[i] = yt2[i];
        }
        h3 = h2;

        real current_value_of_event = event->value(this->t + h3, yt3, dydt3);

        // check if event triggers, if not, return false
        if (current_value_of_event*previous_value_of_event > 0 || previous_value_of_event == 0)
        {
            previous_value_of_event = current_value_of_event;
            return false;
        }
        else if (previous_value_of_event == 0)
        {
            previous_value_of_event = current_value_of_event;
            return true;
        }

        real a = previous_value_of_event, b = current_value_of_event;
        real h_a = 0, h_b = h3;

        for (i = 0; i < MAX_ITERATIONS_SOLVE_EVENT; i++)
        {
            h3 = (h_a*b-h_b*a)/(b-a)*BIAS;
            for (int j = 0; j < this->ode->get_n(); j++)
                yt3[j] = yt[j];
            this->stepper->step_err(t, yt3, h3, err3, dydt, dydt3);
            current_value_of_event = event->value(this->t + h3, yt3, dydt3);
            if (current_value_of_event*previous_value_of_event > 0)
            {
                h_a = h3;
                a = current_value_of_event;
            }
            else
            {
                h_b = h3;
                b = current_value_of_event;
                if( fabs(current_value_of_event) < EVENT_PRECISION)
                {
                    i = MAX_ITERATIONS_SOLVE_EVENT + 1;
                }
            }
        }
       
        if (i == MAX_ITERATIONS_SOLVE_EVENT)
            throw std::runtime_error("precise time of event could not be found");

        previous_value_of_event = current_value_of_event;
        return true;
    }
        
    Integrator::Integrator(ODE &ode, const std::string& stepper_name)
    {
        this->basic_setup();
        this->ode = &ode;
        this->init_stepper(stepper_name);
        this->stepper->set_ODE(ode);
    }

    Integrator::Integrator(ODE &ode, const std::string& stepper_name, const real &atol, const real &rtol)
    {
        this->basic_setup();
        this->ode = &ode;
        this->init_stepper(stepper_name);

        delete this->stepcontroller;
        this->stepcontroller = new StepControllerNR(ode.get_n(), 0, atol, rtol);
    }

    void Integrator::add_event(Event* event)
    {
        switch (event->get_type())
        {
        case EventType::data:
            this->events_data.push_back(event);
            break;

        case EventType::modyfing:
            this->events_modifying.push_back(event);
            break;

        case EventType::terminal:
            this->events_terminal.push_back(event);
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
        for (int i = 0; i < this->ode->get_n(); i++)
            this->yt[i] = this->yt2[i] = y_start[i];
        this->t = t_start;
        this->h = this->h2 = this->h3 = h_start;

        // number of events
        number_of_events_data = events_data.size();
        number_of_events_modifying = events_modifying.size();
        number_of_events_terminal = events_terminal.size();

        // prepare event values
        events_data_values = new real[number_of_events_data];
        events_modifying_values = new real[number_of_events_modifying];
        events_terminal_values = new real[number_of_events_terminal];

        // prepare values of events
        for (int i = 0; i < number_of_events_data; i++)
            events_data_values[i] = events_data[i]->value(t_start, yt, dydt);
        for (int i = 0; i < number_of_events_modifying; i++)
            events_modifying_values[i] = events_modifying[i]->value(t_start, yt, dydt);
        for (int i = 0; i < number_of_events_terminal; i++)
            events_terminal_values[i] = events_terminal[i]->value(t_start, yt, dydt);

        // cycle for calculating new values of y
        while (t < t_end)
        {   
            // std::cout << t << std::endl;
            // taky a step
            this->stepper->step_err(t, yt2, h, err2);

            // null current event
            current_event = nullptr;

            // check modifying events
            for (int i = 0; i < number_of_events_modifying; i++)
            {
                if(this->solve_event(events_modifying[i], events_modifying_values[i]))
                {
                    for (int i = 0; i < this->ode->get_n(); i++)
                    {
                        yt2[i] = yt3[i];
                    }
                    h2 = h3;
                    current_event = events_modifying[i];
                }
            }

            // check terminal events
            for (int i = 0; i < number_of_events_terminal; i++)
            {
                if(this->solve_event(events_terminal[i], events_terminal_values[i]))
                {
                    for (int i = 0; i < this->ode->get_n(); i++)
                    {
                        yt2[i] = yt3[i];
                    }
                    h2 = h3;
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
                    for (int j = 0; j < ode->get_n(); j++)
                        yt2[j] = yt[j];
                    this->stepper->step_err(t, yt2, h2, err2, dydt, dydt2);
                }
            }

            if (i >= MAX_ITERATIONS_HADJUST)
                throw std::runtime_error("optimal step size was not found, MAX_ITERATIONS_HADJUST reached");

            for (int i = 0; i < number_of_events_data; i++)
            {
                if(this->solve_event(events_data[i], events_data_values[i]))
                    events_data[i]->apply(t, yt3, dydt3);
            }

            // "commit" to the step
            for (int i = 0; i < ode->get_n(); i++)
            {
                yt[i] = yt2[i];
                dydt[i] = dydt2[i];
            }
            t += h2;
            if (this->stepcontroller)
                h = h2;
            h2 = h;
            h3 = h;

            // apply event
            if (current_event)
                current_event->apply(t, yt, dydt);

            if (current_event && current_event->get_type() == EventType::terminal)
            {
                break;
            }
        }

        // delete events
        delete[] events_data_values, events_modifying_values, events_terminal_values;
        delete[] yt, yt2, yt3, dydt, dydt2, dydt3, err, err2, err3;
    }
}