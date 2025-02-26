#include <cmath>

#include "gtest/gtest.h"
#include "gravitacek2/integrator/odesystem.hpp"
#include "gravitacek2/integrator/event.hpp"
#include "gravitacek2/integrator/stepperbase.hpp"
#include "gravitacek2/integrator/steppers.hpp"
#include "gravitacek2/integrator/stepcontrollerbase.hpp"
#include "gravitacek2/integrator/stepcontrollers.hpp"

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

class ZeroY : public gr2::Event
{
    protected:
        gr2::real *last_y;
    public:
        ZeroY(gr2::real *last_y) : gr2::Event(gr2::EventType::data){this->last_y = last_y;}
        virtual gr2::real value(const gr2::real &t, const gr2::real y[], const gr2::real dydt[]) override
        {
            return y[1];
        };
        virtual void apply(gr2::real &t, gr2::real y[], gr2::real dydt[])
        {
            this->last_y[0] = y[1];
        };
};

// ========== Event ========== 
TEST(Event, get)
{
    gr2::real y[2] = {1.0, 2.0};
    gr2::real last_y;
    ZeroY event = ZeroY(&last_y);

    EXPECT_EQ(event.get_type(), gr2::EventType::data);
}

TEST(Event, value)
{
    gr2::real y[2] = {1.0, 2.0};
    gr2::real last_y;
    ZeroY event = ZeroY(&last_y);

    EXPECT_EQ(y[1], event.value(0, y, nullptr));
}

TEST(Event, apply)
{
    gr2::real y[2] = {1.0, 2.0};
    gr2::real last_y;
    ZeroY event = ZeroY(&last_y);

    gr2::real t = 0, h = 0;
    event.apply(t, y, nullptr);
    EXPECT_EQ(last_y, y[1]);
}

// ========== ODE ========== 
TEST(ODE, get_n)
{
    // parameters and variables
    gr2::real omega0 = 1.5, xi = 1.0;

    // ODE
    DampedHarmonicOscillator osc = DampedHarmonicOscillator(omega0, xi);

    // test get_n
    ASSERT_EQ(osc.get_n(), 2);
}

TEST(ODE, function)
{
    // parameters and variables
    gr2::real omega0 = 1.5, xi = 1.0;
    gr2::real x0 = 0.5, v0 = 1.5;
    gr2::real y[] = {x0, v0};
    gr2::real dydt[2];
    gr2::real eps = 1e-15;

    // ODE
    DampedHarmonicOscillator osc = DampedHarmonicOscillator(omega0, xi);

    // test function
    osc.function(0, y, dydt);
    ASSERT_NEAR(dydt[0], v0, eps); 
    ASSERT_NEAR(dydt[1], -2*xi*v0-omega0*omega0*x0, eps);
}

TEST(Stepper, DumpedOscillator)
{
    // parameters and variables
    gr2::real omega0 = 2.0, xi = 0.5;
    gr2::real x0 = 1.5, v0 = 0.5;
    gr2::real y[] = {x0, v0};
    gr2::real dydt[2];
    gr2::real h = 0.001;
    gr2::real eps = 1e-7;

    // ODE
    auto osc = std::make_shared<DampedHarmonicOscillator>(omega0, xi);

    // Stepper
    gr2::RK4 stepper = gr2::RK4(false);
    stepper.set_OdeSystem(osc);

    // test order
    ASSERT_EQ(stepper.get_order(), 4);

    // test step
    for (int i = 0; i < int(7/h); i++)
    {
        ASSERT_NEAR(y[0], exactDampedHarmonicOscillator(h*i, omega0, xi, x0, v0), eps);
        stepper.step(i*h, y, h);
    }
}

TEST(Stepper, DumpedOscillatorWithError)
{
    // parameters and variables
    gr2::real omega0 = 2.0, xi = 0.5;
    gr2::real x0 = 1.5, v0 = 0.5;
    gr2::real y[] = {x0, v0};
    gr2::real dydt[2];
    gr2::real err[2];
    gr2::real h = 0.002;
    gr2::real eps = 1e-7;

    // ODE
    auto osc = std::make_shared<DampedHarmonicOscillator>(omega0, xi);

    // Stepper
    gr2::RK4 stepper = gr2::RK4(false);
    stepper.set_OdeSystem(osc);

    // test order
    ASSERT_EQ(stepper.get_order(), 4);

    // test step
    for (int i = 0; i < int(7/h); i++)
    {
        ASSERT_NEAR(y[0], exactDampedHarmonicOscillator(h*i, omega0, xi, x0, v0), eps);
        stepper.step_err(i*h, y, h, err);
    }
}

int main(int argc, char **argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}