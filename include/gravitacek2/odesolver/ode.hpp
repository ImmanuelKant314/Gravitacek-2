#pragma once
#include "gravitacek2/setup.hpp"

namespace gr2
{
    class ODE
    {
    protected:
        int n;
    public:
        ODE(const int &n);
        int get_n() const;
        virtual int function(const gr2::REAL &t, const gr2::REAL y[], gr2::REAL dydt[]) = 0;
    };
}