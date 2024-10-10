#include "gravitacek2/integrator/steppers.hpp"

namespace gr2
{
    RK4::RK4():StepperBase(), k1(nullptr), k2(nullptr), k3(nullptr), k4(nullptr){}

    RK4::~RK4()
    {
        delete[] k1, k2, k3, k4;
    }

    void RK4::set_OdeSystem(OdeSystem& ode)
    {
        int old_n = n;
        this->StepperBase::set_OdeSystem(ode);
        if (old_n != n)
        {
            delete[] k1, k2, k3, k4;
            k1 = new real[n];
            k2 = new real[n];
            k3 = new real[n];
            k4 = new real[n];
        }
    }

    void RK4::reset()
    {
    }

    void RK4::step(const real &t, real y[], const real &h, const real dydt_in[], real dydt_out[]) 
    {
        int i;

        if (dydt_in)
        {
            // copy dydt
            for (i = 0; i < n; i++)
                k1[i] = dydt_in[i];
        }
        else
        {
            // first correction
            for (i = 0; i < n; i++)
                yt[i] = y[i];
            ode->function(t, yt, k1);
        }

        // second correction
        for (i = 0; i < n; i++)
            yt[i] = y[i] + 0.5 * h * k1[i];
        ode->function(t + 0.5 * h, yt, k2);

        // third correction
        for (i = 0; i < n; i++)
            yt[i] = y[i] + 0.5 * h * k2[i];
        ode->function(t + 0.5 * h, yt, k3);

        // fourth correction
        for (i = 0; i < n; i++)
            yt[i] = y[i] + h * k3[i];
        ode->function(t + 0.5 * h, yt, k4);

        // final result
        for (i = 0; i < n; i++)
            y[i] = y[i] + 1.0 / 6 * h * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]);

        if (dydt_out)
        {
            ode->function(t+h, y, dydt_out);
        }
    }

    int RK4::get_order() const
    {
        return 4;
    }

    int RK4::get_err_order() const
    {
        return 5;
    }
} 
