#include "gtest/gtest.h"
#include "gravitacek2/setup.hpp"
#include "gravitacek2/integrator/integrator.hpp"
#include "gravitacek2/integrator/odesystems.hpp"

#include <cmath>
#include <iostream>

gr2::real exactDampedHarmonicOscillator(gr2::real t, gr2::real omega0, gr2::real xi, gr2::real x0, gr2::real v0)
{
    // calculation of omega
    gr2::real omega = sqrtl(omega0*omega0 - xi*xi);

    // calculation of coefficients
    gr2::real A = (v0+xi*x0)/omega;
    gr2::real B = x0;

    // final value
    return expl(-xi*t)*(A*sinl(omega*t) + B*cosl(omega*t));
}

class Bounce : public gr2::Event
{
    protected:
        std::shared_ptr<gr2::DampedHarmonicOscillator> osc;
    public:
        Bounce(std::shared_ptr<gr2::DampedHarmonicOscillator> osc) : gr2::Event(gr2::EventType::modyfing),osc(osc){};
        virtual gr2::real value(const gr2::real &t, const gr2::real &dt, const gr2::real y[], const gr2::real dydt[]) override
        {
            return y[0]-1e-11;
        }
        virtual void apply(gr2::StepperBase* stepper, gr2::real &t, gr2::real &dt, gr2::real y[], gr2::real dydt[]) override
        {
            y[1] = -y[1];
            osc->function(t, y, dydt);
        }
};

class DataMonitoring : public gr2::Event
{
    public:
        std::vector<gr2::real> times;
        std::vector<gr2::real> pos;
        std::vector<gr2::real> vel;
        DataMonitoring() : gr2::Event(gr2::EventType::data)
        {
            times = std::vector<gr2::real>();
            pos = std::vector<gr2::real>();
            vel = std::vector<gr2::real>();
        }
        virtual gr2::real value(const gr2::real &t, const gr2::real &dt, const gr2::real y[], const gr2::real dydt[]) override
        {   
            return 0;
        }
        virtual void apply(gr2::StepperBase* stepper, gr2::real &t, gr2::real &dt, gr2::real y[], gr2::real dydt[])
        {
            times.push_back(t);
            pos.push_back(y[0]);
            vel.push_back(y[1]);
        };
};

class ConstantStepDataMonitoring : public gr2::Event
{
public:
    gr2::real t;
    gr2::real h;
    std::vector<gr2::real> times;
    std::vector<gr2::real> pos;
    std::vector<gr2::real> vel;
    ConstantStepDataMonitoring(gr2::real t_init, gr2::real h) : gr2::Event(gr2::EventType::data)
    {
        t = t_init;
        this->h = h;
        times = std::vector<gr2::real>();
        pos = std::vector<gr2::real>();
        vel = std::vector<gr2::real>();
    }
    virtual gr2::real value(const gr2::real &t, const gr2::real &dt, const gr2::real y[], const gr2::real dydt[]) override
    {   
        return 0;
    }
    virtual void apply(gr2::StepperBase* stepper, gr2::real &t, gr2::real &dt, gr2::real y[], gr2::real dydt[])
    {
        while (this->t<t)
        {
            times.push_back(this->t);
            pos.push_back(stepper->dense_out(0, this->t));
            vel.push_back(stepper->dense_out(1, this->t));
            this->t += h;
        }
    };
};

TEST(Integrator, BouncingDumpedOscilatorNoStepController)
{
    gr2::real omega0 = 1.5, xi = 1.0;
    gr2::real x0 = 0.5, v0 = 1.5;
    gr2::real y0[] = {x0, v0};

    gr2::real eps = 1e-7;
    gr2::real t_start = 0, t_end = 10;
    gr2::real h0 = 0.01;

    auto osc = std::make_shared<gr2::DampedHarmonicOscillator>(omega0, xi);
    auto data = std::make_shared<DataMonitoring>();
    auto bounce = std::make_shared<Bounce>(osc);
    gr2::Integrator integrator = gr2::Integrator(osc, "RK4");

    integrator.add_event(data);
    integrator.add_event(bounce);

    integrator.integrate(y0, t_start, t_end, h0);

    ASSERT_GE(data->times.size(), 10); // check if some data are recorded
    for (int i = 0; i < data->times.size(); i++)
    {
        EXPECT_NEAR(data->pos[i], std::abs(exactDampedHarmonicOscillator(data->times[i], omega0, xi, x0, v0)), eps);
    }
}

TEST(Integrator, BouncingDumberOscilatorStepController)
{
    gr2::real omega0 = 1.5, xi = 1.0;
    gr2::real x0 = 0.5, v0 = 1.5;
    gr2::real y0[] = {x0, v0};

    gr2::real eps = 1e-7;
    gr2::real t_start = 0, t_end = 10;
    gr2::real h0 = 0.01;
    gr2::real atol = 1e-8, rtol = 1e-8;

    auto osc = std::make_shared<gr2::DampedHarmonicOscillator>(omega0, xi);
    auto data = std::make_shared<DataMonitoring>();
    auto bounce = std::make_shared<Bounce>(osc);
    gr2::Integrator integrator = gr2::Integrator(osc, "RK4", atol, rtol);

    integrator.add_event(data);
    integrator.add_event(bounce);

    integrator.integrate(y0, t_start, t_end, h0);

    ASSERT_GE(data->times.size(), 10); // check if some data are recorded
    for (int i = 0; i < data->times.size(); i++)
    {
        EXPECT_NEAR(data->pos[i], abs(exactDampedHarmonicOscillator(data->times[i], omega0, xi, x0, v0)), eps);
    }
}

TEST(Integrator, BouncingDumberOscilatorConstantStepData)
{
    gr2::real omega0 = 1.5, xi = 1.0;
    gr2::real x0 = 0.5, v0 = 1.5;
    gr2::real y0[] = {x0, v0};

    gr2::real eps = 1e-7;
    gr2::real t_start = 0, t_end = 10;
    gr2::real h0 = 0.01;
    gr2::real h_monitor = 0.1;
    gr2::real atol = 1e-8, rtol = 1e-8;

    auto osc = std::make_shared<gr2::DampedHarmonicOscillator>(omega0, xi);
    auto data = std::make_shared<ConstantStepDataMonitoring>(0, h_monitor);
    auto data2 = std::make_shared<DataMonitoring>();
    auto bounce = std::make_shared<Bounce>(osc);
    gr2::Integrator integrator = gr2::Integrator(osc, "RK4", atol, rtol, true);

    integrator.add_event(data);
    integrator.add_event(data2);
    integrator.add_event(bounce);

    integrator.integrate(y0, t_start, t_end, h0);

    ASSERT_GE(data->times.size(), 10); // check if some data are recorded
    for (int i = 0; i < data->times.size(); i++)
    {
        EXPECT_NEAR(data->times[i], h_monitor*i, eps*h_monitor*i);
        EXPECT_NEAR(data->pos[i], std::abs(exactDampedHarmonicOscillator(h_monitor*i, omega0, xi, x0, v0)), eps);
        EXPECT_NEAR(data2->pos[i], std::abs(exactDampedHarmonicOscillator(data2->times[i], omega0, xi, x0, v0)), eps);
    }
}

int main(int argc, char **argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}