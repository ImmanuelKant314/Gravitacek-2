#include <cmath>

#include "gtest/gtest.h"
#include "gravitacek2/integrator/odesystem.hpp"
#include "gravitacek2/integrator/event.hpp"
#include "gravitacek2/integrator/stepperbase.hpp"
#include "gravitacek2/integrator/steppers.hpp"
#include "gravitacek2/integrator/stepcontrollerbase.hpp"
#include "gravitacek2/integrator/stepcontrollers.hpp"

TEST(StepControllerValues, StandardStepController)
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

TEST(StepControllerValues, StepControllerNR)
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

int main(int argc, char **argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}