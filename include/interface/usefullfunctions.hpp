#pragma once

#include "gravitacek2/geomotion/spacetimes.hpp"
#include "gravitacek2/integrator/integrator.hpp"
#include "gravitacek2/integrator/odesystems.hpp"
#include "gravitacek2/chaos/linearized_evolution.hpp"

#include <stdexcept>
#include <iostream>
#include <algorithm>
#include <fstream>
#include <array>
#include <cmath>

class DataRecord : public gr2::Event
{
protected:
    int n;

public:
    std::vector<std::vector<gr2::real>> data;

    DataRecord(int n) : gr2::Event(gr2::EventType::data)
    {
        this->n = n;
    }

    virtual gr2::real value(const gr2::real &t, const gr2::real y[], const gr2::real dydt[]) override
    {
        return 0;
    }

    virtual void apply(gr2::StepperBase* stepper, gr2::real &t, gr2::real y[], gr2::real dydt[]) override
    {
        std::vector<gr2::real> record;
        record.push_back(t);
        for (int i=0; i < n; i++)
            record.push_back(y[i]);
        this->data.push_back(record);
    }
};

class StopOnDisk : public gr2::Event
{
public:
    bool poincare;
    std::vector<std::array<gr2::real, 2>> data;
    StopOnDisk(bool poincare=false) : gr2::Event(gr2::EventType::modyfing), poincare(poincare), data()
    {

    }

    virtual gr2::real value(const gr2::real &t, const gr2::real y[], const gr2::real dydt[]) override
    {
        return y[gr2::Weyl::Z];
    }

    virtual void apply(gr2::StepperBase* stepper, gr2::real &t, gr2::real y[], gr2::real dydt[]) override
    {
        if (poincare)
            data.push_back({y[gr2::Weyl::RHO],y[gr2::Weyl::URHO]});
    }
};

class StopBeforeBlackHole : public gr2::Event
{
public:
    gr2::real rho_min;

    StopBeforeBlackHole(gr2::real rho_min): gr2::Event(gr2::EventType::data, true), rho_min(rho_min)
    {

    }

    virtual gr2::real value(const gr2::real &t, const gr2::real y[], const gr2::real dydt[]) override
    {
        return y[gr2::Weyl::RHO] < rho_min?0:1;
    }

    virtual void apply(gr2::StepperBase* stepper, gr2::real &t, gr2::real y[], gr2::real dydt[]) override
    {
    }
};

class StopTooHighError : public gr2::Event
{
public:
    std::shared_ptr<gr2::Weyl> spt;
    gr2::real E, E_, eps;

    StopTooHighError(std::shared_ptr<gr2::Weyl> spt, gr2::real E, gr2::real eps):gr2::Event(gr2::EventType::data, true), spt(spt), E(E), eps(eps)
    {

    }

    virtual gr2::real value(const gr2::real &t, const gr2::real y[], const gr2::real dydt[]) override
    {
        spt->calculate_metric(y);
        E_ = - spt->get_metric()[gr2::Weyl::T][gr2::Weyl::T]*y[gr2::Weyl::UT];
        return abs(E_-E)/E < eps?1:0;
    }

    virtual void apply(gr2::StepperBase* stepper, gr2::real &t, gr2::real y[], gr2::real dydt[]) override
    {
    }
};