#include "gravitacek2/integrator/stepperbase.hpp"

namespace gr2
{
    StepperBase::StepperBase():t_in(0), h(0), n(0), ode(nullptr), y_in(nullptr), y_out(nullptr), y_err(nullptr), y_help(nullptr), y_cur(nullptr), dydt_in(nullptr), dydt_out(nullptr), dydt_cur(nullptr), dydt_opt(nullptr)
    {}

    StepperBase::~StepperBase()
    {
        delete[] y_in, y_out, y_err, y_cur, y_help;
        delete[] dydt_in, dydt_out, dydt_cur, dydt_opt;
    }

    void StepperBase::set_OdeSystem(std::shared_ptr<OdeSystem> ode)
    {
        int old_n = n;
        n = ode->get_n();
        this->ode = ode;

        // if new ODE has different number of equations
        // create new arrays
        if (old_n != n)
        {
            delete[] y_in, y_out, y_err, y_cur, y_help;
            delete[] dydt_in, dydt_out, dydt_cur, dydt_opt;

            y_in = new real[n];
            y_out = new real[n];
            y_err = new real[n];
            y_cur = new real[n];
            y_help = new real[n];

            dydt_in = new real[n];
            dydt_out = new real[n];
            dydt_cur = new real[n];
            dydt_opt = new real[n];
        }
    }

    void StepperBase::step_err(const real &t, real y[], const real &h, real err[], const bool &dense, const real dydt_in[], real dydt_out[])
    {
        // save time and step internaly
        this->t_in = t;
        this->h = h;

        // copy y to y_in, y_cur and y_err
        for (int i = 0; i < n; i++)
        {
            y_in[i] = y[i];
            y_err[i] = y[i];
            y_cur[i] = y[i];
            y_help[i] = y[i];
        }

        // if dydt_in given, justcopy
        // else calculate it
        if (dydt_in)
            for (int i = 0; i < n; i++)
                this->dydt_in[i] = dydt_in[i];
        else
            ode->function(t, y, this->dydt_in);

        // do two half steps
        this->step(t, y_help, h/2, false, dydt_in=dydt_in, this->dydt_out);
        this->step(t+h/2, y_help, h/2, false, this->dydt_out, dydt_out);

        // do one full step (error estimate)
        this->step(t, y_err, h, false, dydt_in=dydt_in);
        
        // calculate absolute error, copy y_cur to y_out
        for(int i = 0; i < n; i++)
        {
            y_out[i] = y_help[i];
            y_in[i] = y[i];
            err[i] = (y_out[i] - y_err[i]);
            y[i] = y_out[i];
        }

        // if dense, we have to calculate last derivative
        if (dense)
        {
            if (dydt_out)
                for (int i = 0; i < n; i++)
                    this->dydt_out[i] = dydt_out[i];
            else
                this->ode->function(t, y_out, this->dydt_out);
        }
    }

    void StepperBase::prepare_dense()
    {

    }

    real StepperBase::dense_out(const int &i, const real &t)
    {
        real s = (t-t_in)/h;
        real s1 = 1.0-s;
        return s1*y_in[i] + s*y_out[i] - s*s1*((1-2*s)*(y_out[i] - y_in[i]) - s1*h*dydt_in[i] + s*h*dydt_out[i]);
    }
}