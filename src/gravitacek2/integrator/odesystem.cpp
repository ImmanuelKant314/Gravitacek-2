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

    CombinedOdeSystem::CombinedOdeSystem(std::vector<std::shared_ptr<OdeSystem>> odes):OdeSystem(0), odes(odes)
    {
        for (auto &ode:odes)
            this->n+=ode->get_n();
    }

    void CombinedOdeSystem::function(const real &t, const real y[], real dydt[])
    {
        int i = 0;
        for (auto &ode:odes)
        {
            ode->function(t, y+i, dydt+i);
            i+=ode->get_n();
        }
    }
}