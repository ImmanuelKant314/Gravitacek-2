#include "gravitacek2/integrator/odesystem.hpp"

namespace gr2
{
    OdeSystem::OdeSystem(const int &n)
    {
        this->n = n;
    }

    int OdeSystem::get_n() const
    {
        return this->n;
    }
}