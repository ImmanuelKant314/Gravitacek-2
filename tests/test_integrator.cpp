#include "gtest/gtest.h"
#include "gravitacek2/setup.hpp"
#include "gravitacek2/integrator/integrator.hpp"

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

class DampedHarmonicOscillator : public gr2::OdeSystem
{
    protected:
        gr2::real omega0;
        gr2::real xi;

    public:
        DampedHarmonicOscillator(const gr2::real &omega0, const gr2::real &xi):gr2::OdeSystem(2) 
        {
            this->omega0 = omega0;
            this->xi = xi;
        }

        virtual void function(const gr2::real &t, const gr2::real y[], gr2::real dydt[]) override
        {
            dydt[0] = y[1];
            dydt[1] = -2*xi*y[1] - omega0*omega0*y[0];
        }
};

class Bounce : public gr2::Event
{
    protected:
    public:
        Bounce() : gr2::Event(gr2::EventType::modyfing_precise){};
        virtual gr2::real value(const gr2::real &t, const gr2::real y[], const gr2::real dydt[]) override
        {
            return y[0];
        }
        virtual void apply(gr2::real &t, gr2::real y[], gr2::real dydt[]) override
        {
            y[1] = -y[1];
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
        virtual gr2::real value(const gr2::real &t, const gr2::real y[], const gr2::real dydt[]) override
        {
            return 0;
        }
        virtual void apply(gr2::real &t, gr2::real y[], gr2::real dydt[])
        {
            times.push_back(t);
            pos.push_back(y[0]);
            vel.push_back(y[1]);
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

    DataMonitoring data = DataMonitoring();
    Bounce bounce = Bounce();

    DampedHarmonicOscillator osc = DampedHarmonicOscillator(omega0, xi);
    gr2::Integrator integrator = gr2::Integrator(osc, "RK4");

    integrator.add_event(&data);
    integrator.add_event(&bounce);

    integrator.integrate(y0, t_start, t_end, h0);

    ASSERT_GE(data.times.size(), 10); // check if some data are recorded
    for (int i = 0; i < data.times.size(); i++)
    {
        EXPECT_NEAR(data.pos[i], exactDampedHarmonicOscillator(data.times[i], omega0, xi, x0, v0), eps);
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

    DataMonitoring data = DataMonitoring();
    Bounce bounce = Bounce();

    DampedHarmonicOscillator osc = DampedHarmonicOscillator(omega0, xi);
    gr2::Integrator integrator = gr2::Integrator(osc, "RK4", atol, rtol);

    integrator.add_event(&data);
    integrator.add_event(&bounce);

    integrator.integrate(y0, t_start, t_end, h0);

    ASSERT_GE(data.times.size(), 10); // check if some data are recorded
    for (int i = 0; i < data.times.size(); i++)
    {
        EXPECT_NEAR(data.pos[i], exactDampedHarmonicOscillator(data.times[i], omega0, xi, x0, v0), eps);
    }
}

int main(int argc, char **argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}