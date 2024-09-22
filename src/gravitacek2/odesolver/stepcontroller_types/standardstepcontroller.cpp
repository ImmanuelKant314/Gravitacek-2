#include <cmath>

#include "gravitacek2/odesolver/stepcontroller_types.hpp"


namespace gr2
{
    StandardStepController::StandardStepController(const int& n, const int& k, const REAL &eps_abs, const REAL &eps_rel, const REAL &a_y, const REAL &a_dydt, const REAL &S, const REAL &factor):StepController(n)
    {
        this->k = k;
        this->ratios = new REAL[n];
        this->eps_abs = eps_abs;
        this->eps_rel = eps_rel;
        this->a_y = a_y;
        this->a_dydt = a_dydt;
        this->S = S;
        this->factor = factor;
    }

    StandardStepController::~StandardStepController()
    {
        delete[] ratios;
    }

    REAL StandardStepController::hadjust(const REAL y[], const REAL err[], const REAL dydt[], const REAL &h)
    {
        REAL max_ratio = 0;
        REAL h_new = h;
        for (int i = 0; i < n; i++)
        {
            ratios[i] = err[i]/(eps_abs + eps_rel*(a_y*abs(y[i]) + a_dydt*abs(dydt[i])));
            max_ratio = std::max(max_ratio, ratios[i]);
        }

        if (max_ratio > 1.1)
        {
            h_new = h*S*powl(max_ratio, -1.0/k);
            if (h_new*factor < h)
                h_new = h/factor;
        }
        else if(max_ratio < 0.5)
            h_new = h*S*powl(max_ratio, -1.0/(k+1));
            if (h_new > h*factor)
                h_new = h*factor;

        return h_new;
    }
}