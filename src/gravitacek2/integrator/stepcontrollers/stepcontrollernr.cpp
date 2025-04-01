#include <cmath>
// #include <iostream>

#include "gravitacek2/integrator/stepcontrollers.hpp"

namespace gr2
{
    StepControllerNR::StepControllerNR(const int &n, const int &k, const real &atol, const real &rtol, const real &S, const real &factor_decrease, const real &factor_grow) : StepControllerBase(n)
    {
        this->atol = atol;
        this->rtol = rtol;
        this->k = k;
        this->S = S;
        this->factor_decrease = factor_decrease;
        this->factor_grow = factor_grow;
    }

    bool StepControllerNR::hadjust(const real y[], const real err[], const real dydt[], real &h)
    {
        // ========== Calculate error ========== 
        this->err = 0;
        for (int i = 0; i < n; i++)
        {
            scale = atol + rtol*fabs(y[i]);
            // std::cout << "y[i] = " << y[i] << " scale = " << scale << " err[i] = " << err[i] << std::endl;
            real x = (err[i]/scale);
            this->err += x*x;
        }
        // std::cout << "err_ = " << this->err << std::endl;
        this->err = sqrtl(this->err/n);
        // std::cout << "err__ = " << this->err << std::endl;
        // std::cout << "h_old = " << h << std::endl;

        // ========== Calculate step size ==========  
        real h_new = h*S*powl(1.0/this->err, 1.0/this->k);
        if (h_new > factor_grow*h)
            h_new = factor_grow*h;
        else if (h_new < factor_decrease*h)
            h_new = factor_decrease*h;
        h = h_new;
        // std::cout << "h_new = " << h << std::endl;
        
        // ========== Success of step size ========== 
        return this->err <= 1;
    }
}