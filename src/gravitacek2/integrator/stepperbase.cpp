#include "gravitacek2/integrator/stepperbase.hpp"

namespace gr2
{
    StepperBase::StepperBase():t(0), n(0), yt(nullptr), dydt(nullptr), dydt2(nullptr), ode(nullptr)
    {}

    StepperBase::~StepperBase()
    {
        delete[] yt, dydt, dydt2;
    }

    void StepperBase::set_OdeSystem(OdeSystem& ode)
    {
        int old_n = n;
        n = ode.get_n();
        this->ode = &ode;

        // if new ODE has different number of equations
        // create new arrays
        if (old_n != n)
        {
            delete[] yt, dydt, dydt2;
            yt = new real[n];
            yt2 = new real[n];
            dydt = new real[n];
            dydt2 = new real[n];
        }
    }

    void StepperBase::step_err(const real &t, real y[], const real &h, real err[], const real dydt_in[], real dydt_out[])
    {
        // save time internaly
        this->t = t;

        // if dydt_in given, justcopy
        // else calculate it
        if (dydt_in)
            for (int i = 0; i < n; i++)
                dydt[i] = dydt_in[i];
        else
            ode->function(t, y, dydt2);

        // copy y to double calculation
        for (int i = 0; i < n; i++)
            yt2[i] = y[i];

        // do two half steps
        this->step(this->t, y, h/2, dydt_in=dydt2);
        this->step(this->t+h/2, y, h/2);

        // do one full step
        this->step(this->t, yt2, h, dydt_in=dydt2);
        
        // calculate absolute error
        for(int i = 0; i < n; i++)
            err[i] = (y[i] - yt2[i]);

        // if dydt_out given, calculate
        if (dydt_out)
            this->ode->function(this->t+h, y, dydt_out);
    }

}