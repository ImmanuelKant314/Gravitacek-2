// ========== include - my library ========== 
#include "gravitacek2/integrator/integrator.hpp"
#include "gravitacek2/integrator/steppers.hpp"
#include "gravitacek2/integrator/stepcontrollers.hpp"

// ========== include - standard libraries ========== 
#include <stdexcept>
#include <cmath>
// #include <iostream>
// #include <iomanip>

// ========== macros ========== 
#define MAX_ITERATIONS_HADJUST 50 
#define MAX_ITERATIONS_SOLVE_EVENT 30
#define EVENT_PRECISION 1e-8
#define TIME_PRECISION 1e-12

namespace gr2
{
    void Integrator::basic_setup()
    {
        this->ode = nullptr;
        this->stepper = nullptr;
        this->stepcontroller = nullptr;

        this->events_modifying_values = nullptr;

        this->yt = nullptr;
        this->yt2 = nullptr;
        this->yt3 = nullptr;

        this->dydt = nullptr;
        this->dydt2 = nullptr;
        this->dydt3 = nullptr;

        this->err = nullptr;
        this->err2 = nullptr;
        this->err3 = nullptr;

        this->events_data = std::vector<std::shared_ptr<Event>>();
        this->events_modifying = std::vector<std::shared_ptr<Event>>();
    }

    void Integrator::init_stepper(const std::string& stepper_name)
    {
        delete stepper;
        if (stepper_name == "RK4")
            this->stepper = new RK4();
        else if (stepper_name == "DoPr853")
            this->stepper = new DoPr853();
        else
            throw std::invalid_argument("no integrator with given name found");
    }

    bool Integrator::solve_event(std::shared_ptr<Event> event, const real& previous_value_of_event)
    {
        // std::cout << "Solve event" << std::endl;
        // prepare values
        int i;

        // calculate current value of event
        real current_value_of_event = event->value(t2, h2, yt2, dydt2);
        // std::cout << previous_value_of_event << " " << current_value_of_event << std::endl;

        // check if event is triggered
        // if event is not triggered return false 
        // else if current value of event is close enough, return true
        if (current_value_of_event*previous_value_of_event > 0 || previous_value_of_event == 0)
            return false;
        // else if (fabs(previous_value_of_event) < EVENT_PRECISION) 
        // {
        //     std::cout << "Previous value of event too small" << std::endl;
        //     return false;
        // }

        // prepare values for secant method
        real a = previous_value_of_event, b = current_value_of_event;
        real h_a = 0, h_b = h2;

        // std::cout << "==========" << std::endl;
        // std::cout << current_value_of_event << std::endl;
        // std::cout << h2 << std::endl;

        // ========== Secant method ==========
        // std::cout << "Let us iterate" << std::endl;
        for (i = 0; i < MAX_ITERATIONS_SOLVE_EVENT; i++)
        {
            gr2::real avg = 0.5*(h_a+h_b);
            h3 = avg + 0.8*((h_a*b-h_b*a)/(b-a)-avg); // new value of step size

            // copy starting value of yt to yt3
            for (int j = 0; j < this->ode->get_n(); j++)
                yt3[j] = yt[j];

            // take step and calculate new value of event
            this->stepper->step_err(t, yt3, h3, err3, dense, dydt, dydt3);
            t3 = t + h3;
            current_value_of_event = event->value(t3, h3, yt3, dydt3);
            // std::cout << std::scientific;
            // std::cout << "current value of event = " << current_value_of_event << std::endl;

            // reduce interval for h
            if (current_value_of_event*a > 0)
            {
                h_a = h3;
                a = current_value_of_event;
            }
            else
            {
                h_b = h3;
                b = current_value_of_event;
            }
            // if new value of event is close enought, stop for-cycle
            if( (h_b-h_a) < TIME_PRECISION*std::max(h_a, h_b))
            {
                // std::cout << "precision of time"<< std::endl;
                return true;
            }
            if (std::abs(current_value_of_event) < EVENT_PRECISION )
            {
                // std::cout << "precision of event" << std::endl;
                return true;
            }
        }
        if (i == MAX_ITERATIONS_SOLVE_EVENT)
        {
            // std::cout << std::scientific;
            // for (int ii = 0; ii < 8; ii++)
            //     std::cout << this->yt[ii] << " ";
            // for (int ii = 0; ii < 8; ii++)
            //     std::cout << this->yt3[ii] << " ";
            // std::cout << std::endl;
            // std::cout << current_value_of_event << std::endl;
            // std::cout << h2 << std::endl;
            // std::cout << h_a << " " << h_b << " " << h_b-h_a << std::endl;
            // std::cout << a << " " << b << " " << b-a << std::endl;
            // std::cin.get();
            // std::cout << "z = " << yt[3] << std::endl;
            // exit(1);
            throw std::runtime_error("precise time of event could not be found");
        }
        return true;
    }
        
    Integrator::Integrator(std::shared_ptr<OdeSystem> ode, const std::string& stepper_name, const bool &dense)
    {
        this->basic_setup();
        this->dense = dense;
        this->ode = ode;
        this->init_stepper(stepper_name);
        this->stepper->set_OdeSystem(ode);

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
    }

    Integrator::Integrator(std::shared_ptr<OdeSystem> ode, const std::string& stepper_name, const real &atol, const real &rtol, const bool &dense)
    {
        this->basic_setup();
        this->ode = ode;
        this->dense = dense;
        this->init_stepper(stepper_name);
        this->stepper->set_OdeSystem(ode);
        this->stepcontroller = new StepControllerNR(ode->get_n(), this->stepper->get_err_order(), atol, rtol, 0.8, 0.2, 10.0);
        
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
        this->events_modifying_values = new real[n];
    }

    // Integrator::Integrator(std::shared_ptr<OdeSystem> ode, const std::string& stepper_name, const real &atol, const real &rtol, const real &min_step, const real &max_step, const bool &dense = false): h_boundary(true)
    // {
    //     this->basic_setup();
    //     this->ode = ode;
    //     this->dense = dense;
    //     this->init_stepper(stepper_name);
    //     this->stepper->set_OdeSystem(ode);
    //     delete this->stepcontroller;
    //     this->stepcontroller = new StepControllerNR(ode->get_n(), this->stepper->get_err_order(), atol, rtol, 0.95, 0.2, 10.0);
    // }

    void Integrator::add_event(std::shared_ptr<Event> event)
    {
        switch (event->get_type())
        {
        case EventType::data:
            this->events_data.push_back(event);
            break;
        case EventType::modyfing:
            this->events_modifying.push_back(event);
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
        delete[] yt;
        delete[] yt2;
        delete[] yt3;
        delete[] dydt;
        delete[] dydt2;
        delete[] dydt3;
        delete[] err;
        delete[] err2;
        delete[] err3;
    }

    void Integrator::integrate(const real y_start[], const real &t_start, const real &t_end, const real &h_start)
    {
        // prepare variables
        int i;
        int n = this->ode->get_n();


        // copy values internaly
        for (int i = 0; i < n; i++)
            this->yt[i] = this->yt2[i] = y_start[i];
        this->t = t_start;
        this->h = this->h2 = this->h3 = h_start;
        t2 = t3 = t;
        this->ode->function(t, yt, dydt);

        // number of events
        number_of_events_modifying = events_modifying.size();

        // prepare event values
        events_modifying_values = new real[number_of_events_modifying];

        // prepare values of events
        for (int i = 0; i < number_of_events_modifying; i++)
            events_modifying_values[i] = events_modifying[i]->value(t_start, h, yt, dydt);

        // cycle for calculating new values of y
        while (t < t_end)
        {   
            // taky a step
            this->stepper->step_err(t, yt2, h, err2, dense, dydt, dydt2);
            t2 = t + h;

            // null current event
            current_event = nullptr;
            current_event_terminal = false;

            // check modifying events
            int current_event_index = 0;
            for (int i = 0; i < number_of_events_modifying; i++)
            {
                if(this->solve_event(events_modifying[i], events_modifying_values[i]))
                {
                    for (int i = 0; i < n; i++)
                    {
                        yt2[i] = yt3[i];
                        dydt2[i] = dydt3[i];
                        err2[i] = err3[i];
                    }
                    h2 = h3;
                    t2 = t3;
                    current_event = events_modifying[i];
                    current_event_terminal = current_event->get_terminal();
                    current_event_index = i;
                    // std::cout << "Event activated" << std::endl;
                }
            }

            i = 0;
            if (this->stepcontroller)
            {
                // std::cout << std::scientific;
                // std::cout << "Lets do it" << std::endl;
                // std::cout << "z = " << yt[3] << std::endl;
                // std::cout << "z = " << yt2[3] << std::endl;
                // std::cout << "err = " << err2[3] << std::endl;
                // std::cout << "h2 = " << h2 << std::endl;
                for (i = 0; i < MAX_ITERATIONS_HADJUST; i++)
                {
                    if(this->stepcontroller->hadjust(this->yt2, this->err2, this->dydt2, this->h2))
                        break;
                    for (int j = 0; j < n; j++)
                        yt2[j] = yt[j];
                    this->stepper->step_err(t, yt2, h2, err2, dense, dydt, dydt2);
                    t2 = t + h2;

                    // check current event
                    // std::cout << current_event << std::endl;
                    if (current_event && (events_modifying_values[current_event_index]*current_event->value(t2, h2, yt2, dydt2) <= 0))
                    {
                        // std::cout << t2 << " not delete event " << yt2[3] << std::endl;
                        // std::cout << current_event->value(t2, yt2, dydt2) << std::endl;
                    }
                    else
                    {
                        // std::cout << t2 << " delete event " << yt2[3] << std::endl;
                        // std::cout << events_modifying_values[0] << " " << events_modifying[0]->value(t2, yt2, dydt2) << std::endl;
                        // if (std::abs(yt2[3]) < 1e-5)
                        //     throw;
                        current_event = nullptr;
                        current_event_terminal = false;
                    }

                    // current_event = nullptr;
                    // current_event_terminal = false;
                }
            }

            // if too much iteration for hadjust, throw exception
            if (i >= MAX_ITERATIONS_HADJUST)
            {
                // std::cout << h2 << std::endl;
                // std::cout << yt[3] << std::endl;
                // exit(1);
                throw std::runtime_error("optimal step size was not found, MAX_ITERATIONS_HADJUST reached");
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
            // t2 = t3 = t + h;

            if (dense)
                this->stepper->prepare_dense();

            // apply event
            if (current_event)
                current_event->apply(stepper, t, h, yt, dydt);
            
            // data events
            for (auto &event : events_data)
                if (event->value(t, h, yt, dydt) == 0)
                {
                    event->apply(stepper, t, h, yt, dydt);
                    current_event_terminal = std::max(event->get_terminal(), current_event_terminal);
                }

            // kill precise event
            if (current_event_terminal)
                break;

            // calculate new values of events
            for (int i = 0; i < number_of_events_modifying; i++)
            {
                events_modifying_values[i] = events_modifying[i]->value(t, h, yt, dydt);
            }

            // propagate changes from events
            for (int i = 0; i < n; i++)
            {
                yt2[i] = yt[i];
                dydt2[i] = dydt[i];
            }
        }
        delete[] events_modifying_values;
    }
}