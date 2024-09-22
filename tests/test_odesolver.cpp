#include <cmath>

#include "gtest/gtest.h"
#include "gravitacek2/odesolver/ode.hpp"
#include "gravitacek2/odesolver/stepper.hpp"
#include "gravitacek2/odesolver/steper_types.hpp"
#include "gravitacek2/odesolver/stepcontroller.hpp"
#include "gravitacek2/odesolver/stepcontroller_types.hpp"

gr2::REAL exactDampedHarmonicOscillator(gr2::REAL t, gr2::REAL omega0, gr2::REAL xi, gr2::REAL x0, gr2::REAL v0)
{
    // calculation of omega
    gr2::REAL omega = sqrtl(omega0*omega0 - xi*xi);

    // calculation of coefficients
    gr2::REAL A = (v0+xi*x0)/omega;
    gr2::REAL B = x0;

    // final value
    return expl(-xi*t)*(A*sinl(omega*t) + B*cosl(omega*t));
}

class DampedHarmonicOscillator : public gr2::ODE
{
    protected:
        gr2::REAL omega0;
        gr2::REAL xi;

    public:
        DampedHarmonicOscillator(const gr2::REAL &omega0, const gr2::REAL &xi):gr2::ODE(2) 
        {
            this->omega0 = omega0;
            this->xi = xi;
        }

        virtual void function(const gr2::REAL &t, const gr2::REAL y[], gr2::REAL dydt[]) override
        {
            dydt[0] = y[1];
            dydt[1] = -2*xi*y[1] - omega0*omega0*y[0];
        }
};

TEST(odesolver_test, standardtestcontroller)
{
    int n = 2, k = 4;
    gr2::REAL eps_abs = 1e-10, eps_rel = 2e-10, a_y = 1, a_dydt = 1;
    gr2::REAL h_old = 1e-3, h_new;

    gr2::REAL y[2] = {0.5, 0.5};
    gr2::REAL dydt[2] = {0.1, 0.1};
    gr2::REAL err[2];

    gr2::REAL eps = 1e-5;

    gr2::StandardStepController stepcontroller = gr2::StandardStepController(n, k, eps_abs, eps_rel, a_y, a_dydt);

    // large error
    err[0] = err[1] = 1e-3;
    h_new = stepcontroller.hadjust(y, err, dydt, h_old);
    EXPECT_NEAR(h_new/h_old, 1.0/5, eps);

    // little bit larger error
    err[0] = err[1] = 1.2e-10;
    h_new = stepcontroller.hadjust(y, err, dydt, h_old);
    EXPECT_NEAR(h_new/h_old, 0.907671, eps);

    // normal error 
    err[0] = err[1] = 1e-10;
    h_new = stepcontroller.hadjust(y, err, dydt, h_old);
    EXPECT_NEAR(h_new/h_old, 1.0, eps);

    // little bit small error
    err[0] = err[1] = 0.4e-10;
    h_new = stepcontroller.hadjust(y, err, dydt, h_old);
    EXPECT_NEAR(h_new/h_old, 1.14107, eps);

    // small error
    err[0] = err[1] = 1e-15;
    h_new = stepcontroller.hadjust(y, err, dydt, h_old);
    EXPECT_NEAR(h_new/h_old, 5, eps);
}

TEST(odesolver_test, damped_harmonic_oscillator_ode)
{
    // parameters and variables
    gr2::REAL omega0 = 1.5, xi = 1.0;
    gr2::REAL x0 = 0.5, v0 = 1.5;
    gr2::REAL y[] = {x0, v0};
    gr2::REAL dydt[2];
    gr2::REAL eps = 1e-15;

    // ODE
    DampedHarmonicOscillator osc = DampedHarmonicOscillator(omega0, xi);

    // test function
    osc.function(0, y, dydt);
    ASSERT_NEAR(dydt[0], v0, eps); 
    ASSERT_NEAR(dydt[1], -2*xi*v0-omega0*omega0*x0, eps);

    // test get_n
    ASSERT_EQ(osc.get_n(), 2);
}

TEST(odesolver_test, damped_harmonic_oscillator_stepper)
{
    // parameters and variables
    gr2::REAL omega0 = 2.0, xi = 0.5;
    gr2::REAL x0 = 1.5, v0 = 0.5;
    gr2::REAL y[] = {x0, v0};
    gr2::REAL dydt[2];
    gr2::REAL h = 0.001;
    gr2::REAL eps = 1e-7;

    // ODE
    DampedHarmonicOscillator osc = DampedHarmonicOscillator(omega0, xi);

    // Stepper
    gr2::RK4 stepper = gr2::RK4();
    stepper.set_ODE(osc);

    // test order
    ASSERT_EQ(stepper.get_order(), 4);

    // test step
    for (int i = 0; i < int(7/h); i++)
    {
        ASSERT_NEAR(y[0], exactDampedHarmonicOscillator(h*i, omega0, xi, x0, v0), eps);
        stepper.step(i*h, y, h);
    }
}

TEST(odesolver_test, damped_harmonic_oscillator_stepper_with_err)
{
    // parameters and variables
    gr2::REAL omega0 = 2.0, xi = 0.5;
    gr2::REAL x0 = 1.5, v0 = 0.5;
    gr2::REAL y[] = {x0, v0};
    gr2::REAL dydt[2];
    gr2::REAL err[2];
    gr2::REAL h = 0.002;
    gr2::REAL eps = 1e-7;

    // ODE
    DampedHarmonicOscillator osc = DampedHarmonicOscillator(omega0, xi);

    // Stepper
    gr2::RK4 stepper = gr2::RK4();
    stepper.set_ODE(osc);

    // test order
    ASSERT_EQ(stepper.get_order(), 4);

    // test step
    for (int i = 0; i < int(7/h); i++)
    {
        ASSERT_NEAR(y[0], exactDampedHarmonicOscillator(h*i, omega0, xi, x0, v0), eps);
        stepper.step(i*h, y, h, err);
    }
}
int main(int argc, char **argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}