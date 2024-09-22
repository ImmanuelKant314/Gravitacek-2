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
        virtual REAL hadjust(const REAL y[], const REAL err[], const REAL dydt[], const REAL &h) = 0;
    };
}