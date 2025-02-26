#include "gravitacek2/integrator/steppers.hpp"

namespace gr2
{
    RK4::RK4():StepperBase(), k1(nullptr), k2(nullptr), k3(nullptr), k4(nullptr){}

    RK4::~RK4()
    {
        delete[] k2, k3, k4;
    }

    void RK4::set_OdeSystem(std::shared_ptr<OdeSystem> ode)
    {
        int old_n = n;
        this->StepperBase::set_OdeSystem(ode);
        if (old_n != n)
        {
            delete[] k2, k3, k4;
            k1 = dydt_in;
            k2 = new real[n];
            k3 = new real[n];
            k4 = new real[n];
        }
    }

    void RK4::reset()
    {
    }

    void RK4::step(const real &t, real y[], const real &h, const bool &dense, const real dydt_in[], real dydt_out[]) 
    {
        int i;
        
        // save time and step internaly
        this->t_in = t;
        this->h = h;

        // copy y to y_in, y_cur
        for (int i = 0; i < n; i++)
        {
            y_in[i] = y[i];
            y_cur[i] = y[i];
        }

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
                y_cur[i] = y_in[i];
            ode->function(t, y_cur, k1);
        }

        // second correction
        for (i = 0; i < n; i++)
            y_cur[i] = y_in[i] + 0.5 * h * k1[i];
        ode->function(t + 0.5 * h, y_cur, k2);

        // third correction
        for (i = 0; i < n; i++)
            y_cur[i] = y_in[i] + 0.5 * h * k2[i];
        ode->function(t + 0.5 * h, y_cur, k3);

        // fourth correction
        for (i = 0; i < n; i++)
            y_cur[i] = y_in[i] + h * k3[i];
        ode->function(t + 0.5 * h, y_cur, k4);

        // final result
        for (i = 0; i < n; i++)
        {
            y_out[i] = y_in[i] + 1.0 / 6 * h * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]);
            y[i] = y_out[i];
        }

        if (dydt_out)
        {
            ode->function(t+h, y, dydt_out);
            if (dense)
                for (i = 0; i < n; i++)
                    this->dydt_out[i] = dydt_out[i];
        }
        else if (dense)
        {
            ode->function(t+h, y, this->dydt_out);
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
