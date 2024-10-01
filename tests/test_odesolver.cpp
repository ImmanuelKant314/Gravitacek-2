#include <cmath>

#include "gtest/gtest.h"
#include "gravitacek2/odesolver/ode.hpp"
#include "gravitacek2/odesolver/event.hpp"
#include "gravitacek2/odesolver/stepper.hpp"
#include "gravitacek2/odesolver/steper_types.hpp"
#include "gravitacek2/odesolver/stepcontroller.hpp"
#include "gravitacek2/odesolver/stepcontroller_types.hpp"

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

class DampedHarmonicOscillator : public gr2::ODE
{
    protected:
        gr2::real omega0;
        gr2::real xi;

    public:
        DampedHarmonicOscillator(const gr2::real &omega0, const gr2::real &xi):gr2::ODE(2) 
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
        virtual gr2::real value(const gr2::real &t, const gr2::real y[], const gr2::real err[], const gr2::real dydt[], const gr2::real &h) override
        {
            return y[1];
        };
        virtual void apply(gr2::real &t, gr2::real y[], gr2::real err[], gr2::real dydt[], gr2::real &h)
        {
            this->last_y[0] = y[1];
        };
};

TEST(step_controller, standard_step_controller)
{
    int n = 2, k = 4;
    gr2::real eps_abs = 1e-10, eps_rel = 2e-10, a_y = 1, a_dydt = 1;
    gr2::real h_old = 1e-4, h_new;
    bool test;

    gr2::real y[2] = {0.5, 0.5};
    gr2::real dydt[2] = {0.1, 0.1};
    gr2::real err[2];

    gr2::real eps = 1e-3;

    gr2::StandardStepController stepcontroller = gr2::StandardStepController(n, k, eps_abs, eps_rel, a_y, a_dydt);

    // large error
    err[0] = err[1] = 1e-4;
    h_new = h_old; 
    test = stepcontroller.hadjust(y, err, dydt, h_new);
    EXPECT_FALSE(test);
    EXPECT_NEAR(h_new/h_old, 1.0/5, eps);

    // little bit larger error
    err[0] = err[1] = 1.2e-10;
    h_new = h_old;
    test = stepcontroller.hadjust(y, err, dydt, h_new);
    EXPECT_FALSE(test);
    EXPECT_NEAR(h_new/h_old, 0.907671, eps);

    // normal error 
    err[0] = err[1] = 1e-10;
    h_new = h_old;
    test = stepcontroller.hadjust(y, err, dydt, h_new);
    EXPECT_TRUE(test);
    EXPECT_NEAR(h_new/h_old, 1.0, eps);

    // little bit small error
    err[0] = err[1] = 0.4e-10;
    h_new = h_old;
    test = stepcontroller.hadjust(y, err, dydt, h_new);
    EXPECT_TRUE(test);
    EXPECT_NEAR(h_new/h_old, 1.14107, eps);

    // small error
    err[0] = err[1] = 1e-15;
    h_new = h_old;
    test = stepcontroller.hadjust(y, err, dydt, h_new);
    EXPECT_TRUE(test);
    EXPECT_NEAR(h_new/h_old, 5, eps);
}

TEST(step_controller, step_controller_nr)
{
    // ========== Parameters for controller ========== 
    int n = 2, k = 4;
    gr2::real atol = 1e-10, rtol = 2e-10;
    gr2::real h_old = 1e-4, h_new;

    gr2::StepControllerNR stepcontroller = gr2::StepControllerNR(n, k, atol, rtol);

    // ========== Test example ========== 
    gr2::real y[2] = {0.5, 0.5};
    gr2::real dydt[2] = {0.1, 0.1};
    gr2::real err[2];

    // ========== Variables for test ========== 
    gr2::real eps = 1e-3;
    bool test;

    // ========== Testing ========== 
    // big error
    err[0] = err[1] = 1e-3;
    h_new = h_old; 
    test = stepcontroller.hadjust(y, err, dydt, h_new);
    EXPECT_FALSE(test);
    EXPECT_NEAR(h_new/h_old, 1.0/5, eps);

    // little big error
    err[0] = err[1] = 1e-9;
    h_new = h_old; 
    test = stepcontroller.hadjust(y, err, dydt, h_new);
    EXPECT_FALSE(test);
    EXPECT_NEAR(h_new/h_old, 0.63530, eps);

    // little small error
    err[0] = err[1] = 1e-11;
    h_new = h_old; 
    test = stepcontroller.hadjust(y, err, dydt, h_new);
    EXPECT_TRUE(test);
    EXPECT_NEAR(h_new/h_old, 2.0090, eps);

    // small error
    err[0] = err[1] = 1e-15;
    h_new = h_old; 
    test = stepcontroller.hadjust(y, err, dydt, h_new);
    EXPECT_TRUE(test);
    EXPECT_NEAR(h_new/h_old, 10, eps);
}

TEST(events, zero_y)
{
    gr2::real y[2] = {1.0, 2.0};
    gr2::real last_y;
    ZeroY event = ZeroY(&last_y);

    EXPECT_EQ(event.get_type(), gr2::EventType::data);
    EXPECT_EQ(y[1], event.value(0, y, nullptr, nullptr, 0));
    gr2::real t = 0, h = 0;
    event.apply(t, y, nullptr, nullptr, h);
    EXPECT_EQ(last_y, y[1]);
}

TEST(odesolver_test, damped_harmonic_oscillator_ode)
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

    // test get_n
    ASSERT_EQ(osc.get_n(), 2);
}

TEST(odesolver_test, damped_harmonic_oscillator_stepper)
{
    // parameters and variables
    gr2::real omega0 = 2.0, xi = 0.5;
    gr2::real x0 = 1.5, v0 = 0.5;
    gr2::real y[] = {x0, v0};
    gr2::real dydt[2];
    gr2::real h = 0.001;
    gr2::real eps = 1e-7;

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
    gr2::real omega0 = 2.0, xi = 0.5;
    gr2::real x0 = 1.5, v0 = 0.5;
    gr2::real y[] = {x0, v0};
    gr2::real dydt[2];
    gr2::real err[2];
    gr2::real h = 0.002;
    gr2::real eps = 1e-7;

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
        stepper.step_err(i*h, y, h, err);
        std::cout << err[0] << ", " << err[1] << std::endl;
    }
}
int main(int argc, char **argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}