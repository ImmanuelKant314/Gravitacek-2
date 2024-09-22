#pragma once
#include "gravitacek2/setup.hpp"

namespace gr2
{
    class StepController
    {
    protected:
        int n;
    public:
        StepController(const int &n);
        virtual real hadjust(const real y[], const real err[], const real dydt[], const real &h) = 0;
    };
}