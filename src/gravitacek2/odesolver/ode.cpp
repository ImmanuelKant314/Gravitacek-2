#include "gravitacek2/odesolver/ode.hpp"

namespace gr2
{
    ODE::ODE(const int &n)
    {
        this->n = n;
    }

    int ODE::get_n() const
    {
        return this->n;
    }
}