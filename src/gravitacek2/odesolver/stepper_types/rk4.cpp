#include "gravitacek2/odesolver/steper_types.hpp"

namespace gr2
{
    RK4::RK4():Stepper(), k1(nullptr), k2(nullptr), k3(nullptr), k4(nullptr){}

    RK4::~RK4()
    {
        delete[] k1, k2, k3, k4;
    }

    void RK4::set_ODE(ODE& ode)
    {
        int old_n = n;
        this->Stepper::set_ODE(ode);
        if (old_n != n)
        {
            delete[] k1, k2, k3, k4;
            k1 = new REAL[n];
            k2 = new REAL[n];
            k3 = new REAL[n];
            k4 = new REAL[n];
        }
    }

    void RK4::reset()
    {
    }

    void RK4::step(const REAL &t, REAL y[], const REAL &h, const REAL dydt_in[], REAL dydt_out[]) 
    {
        int i;

        // first correction
        for (i = 0; i < n; i++)
            yt[i] = y[i];
        ode->function(t, yt, k1);

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
    }

    int RK4::get_order() const
    {
        return 4;
    }
} 
