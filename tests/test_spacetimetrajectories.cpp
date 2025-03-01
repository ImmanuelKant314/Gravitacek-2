#include <iostream>
#include <fstream>
#include <cmath>

#include "gtest/gtest.h"
#include "gravitacek2/setup.hpp"
#include "gravitacek2/integrator/steppers.hpp"
#include "gravitacek2/geomotion/geomotion.hpp"
#include "gravitacek2/geomotion/spacetimes.hpp"

TEST(SchwarzschildTrajectory, IntegralsOfMotionPlanar)
{
    gr2::real eps = 1e-13;

    // prepare objects
    auto stepper = std::make_shared<gr2::DoPr853>();
    auto spt = std::make_shared<gr2::Schwarzschild>(1.0);
    stepper->set_OdeSystem(spt);
    gr2::real y[8]{};

    // initial conditions - position
    y[gr2::Schwarzschild::R] = 16;
    y[gr2::Schwarzschild::THETA] = gr2::pi_2;

    // initial conditions - velocity
    spt->calculate_metric(y);
    gr2::real L = 3.6823981191047921;
    y[gr2::Schwarzschild::UPHI] = L/spt->get_metric()[gr2::Schwarzschild::PHI][gr2::Schwarzschild::PHI];
    y[gr2::Schwarzschild::UT] = sqrtl((-1 - L*y[gr2::Schwarzschild::UPHI])/spt->get_metric()[gr2::Schwarzschild::T][gr2::Schwarzschild::T]);
    gr2::real E = -spt->get_metric()[gr2::Schwarzschild::T][gr2::Schwarzschild::T]*y[gr2::Schwarzschild::UT];
    ASSERT_NEAR(E, 0.9598686615055122147, 1e-15); // check if energy is correct

    // initial conditions - steps
    gr2::real dt = 0.2;
    int N = 5000;
    gr2::real r_min = 16;

    for (int i = 0; i < N; i++)
    {
        stepper->step(0, y, dt, false, nullptr, nullptr);
        spt->calculate_metric(y);
        // energy
        EXPECT_NEAR(E, -spt->get_metric()[gr2::Schwarzschild::T][gr2::Schwarzschild::T]*y[gr2::Schwarzschild::UT], eps);
        // angular momentum
        EXPECT_NEAR(L, spt->get_metric()[gr2::Schwarzschild::PHI][gr2::Schwarzschild::PHI]*y[gr2::Schwarzschild::UPHI], eps);// norm of four-velocity
        gr2::real norm2 = 0;
        for (int j = 0; j < 4; j++)
            norm2 += spt->get_metric()[j][j]*y[4+j]*y[4+j];
        EXPECT_NEAR(-1, norm2, eps);
        r_min = std::min(r_min, y[gr2::Schwarzschild::R]);
    }
    std::cout << r_min << std::endl;
    std::cout <<spt->get_metric()[gr2::Schwarzschild::PHI][gr2::Schwarzschild::PHI]*y[gr2::Schwarzschild::UPHI] - L << std::endl;
}

TEST(SchwarzschildTrajectory, IntegralsOfMotionGeneral)
{
    gr2::real eps = 1e-13;

    // prepare objects
    auto stepper = std::make_shared<gr2::DoPr853>();
    auto spt = std::make_shared<gr2::Schwarzschild>(1.0);
    stepper->set_OdeSystem(spt);
    gr2::real y[8]{};

    // initial conditions - position
    y[gr2::Schwarzschild::R] = 16;
    y[gr2::Schwarzschild::THETA] = gr2::pi_2;

    // initial conditions - velocity
    spt->calculate_metric(y);
    gr2::real L = 3.6823981191047921;
    gr2::real E = 0.97;
    y[gr2::Schwarzschild::UPHI] = L/spt->get_metric()[gr2::Schwarzschild::PHI][gr2::Schwarzschild::PHI];
    y[gr2::Schwarzschild::UT] = -E/spt->get_metric()[gr2::Schwarzschild::T][gr2::Schwarzschild::T];
    y[gr2::Schwarzschild::UR] = 0.0;

    gr2::real theta = 0;

    gr2::real norm2 = 0;
    for (int j = 0; j < 4; j++)
        norm2 += spt->get_metric()[j][j]*y[4+j]*y[4+j];
    ASSERT_GT(-1-norm2, 0);
    y[gr2::Schwarzschild::UTHETA] = sqrtl((-1-norm2)/spt->get_metric()[gr2::Schwarzschild::THETA][gr2::Schwarzschild::THETA]);

    // initial conditions - steps
    const gr2::real dt = 0.7;
    int N = 10000;

    for (int i = 0; i < N; i++)
    {
        stepper->step(0, y, dt);
        spt->calculate_metric(y);
        // energy
        EXPECT_NEAR(E, -spt->get_metric()[gr2::Schwarzschild::T][gr2::Schwarzschild::T]*y[gr2::Schwarzschild::UT], eps);
        // angular momentum
        EXPECT_NEAR(L, spt->get_metric()[gr2::Schwarzschild::PHI][gr2::Schwarzschild::PHI]*y[gr2::Schwarzschild::UPHI], eps);
        // norm
        gr2::real norm2 = 0;
        for (int j = 0; j < 4; j++)
            norm2 += spt->get_metric()[j][j]*y[4+j]*y[4+j];
        EXPECT_NEAR(-1, norm2, eps);

        theta = std::max(theta, y[gr2::Schwarzschild::THETA]);
    }
}

TEST(WeylSchwarzschildTrajectory, IntegralsOfMotionPlanar)
{
    gr2::real eps = 1e-13;

    // prepare objects
    auto stepper = std::make_shared<gr2::DoPr853>();
    auto spt = std::make_shared<gr2::WeylSchwarzschild>(1.0, gr2::exact, gr2::exact);
    stepper->set_OdeSystem(spt);
    gr2::real y[8]{};

    // initial conditions - position
    y[gr2::Weyl::RHO] = sqrtl(16*14);

    // initial conditions - velocity
    spt->calculate_metric(y);
    gr2::real L = 3.6823981191047921;
    y[gr2::Weyl::UPHI] = L/spt->get_metric()[gr2::Weyl::PHI][gr2::Weyl::PHI];
    y[gr2::Weyl::UT] = sqrtl((-1 - L*y[gr2::Weyl::UPHI])/spt->get_metric()[gr2::Weyl::T][gr2::Weyl::T]);
    gr2::real E = -spt->get_metric()[gr2::Weyl::T][gr2::Weyl::T]*y[gr2::Weyl::UT];
    ASSERT_NEAR(E, 0.9598686615055122147, 1e-15); // check if energy is correct

    // initial conditions - steps
    gr2::real dt = 0.2;
    int N = 5000;
    gr2::real r_min = 16;

    for (int i = 0; i < N; i++)
    {
        stepper->step(0, y, dt, false, nullptr, nullptr);
        spt->calculate_metric(y);
        // energy
        EXPECT_NEAR(E, -spt->get_metric()[gr2::Weyl::T][gr2::Weyl::T]*y[gr2::Weyl::UT], eps);
        // angular momentum
        EXPECT_NEAR(L, spt->get_metric()[gr2::Weyl::PHI][gr2::Weyl::PHI]*y[gr2::Weyl::UPHI], eps);// norm of four-velocity
        gr2::real norm2 = 0;
        for (int j = 0; j < 4; j++)
            norm2 += spt->get_metric()[j][j]*y[4+j]*y[4+j];
        EXPECT_NEAR(-1, norm2, eps);
        r_min = std::min(r_min, y[gr2::Weyl::RHO]);
    }
    std::cout << r_min << std::endl;
}

TEST(WeylSchwarzschildTrajectory, IntegralsOfMotionGeneral)
{
    gr2::real eps = 1e-13;

    // prepare objects
    auto stepper = std::make_shared<gr2::DoPr853>();
    auto spt = std::make_shared<gr2::WeylSchwarzschild>(1.0, gr2::exact, gr2::exact);
    stepper->set_OdeSystem(spt);
    gr2::real y[8]{};

    // initial conditions - position
    y[gr2::Weyl::RHO] = sqrtl(16*14);

    // initial conditions - velocity
    spt->calculate_metric(y);
    gr2::real L = 3.6823981191047921;
    gr2::real E = 0.97;
    y[gr2::Weyl::UPHI] = L/spt->get_metric()[gr2::Weyl::PHI][gr2::Weyl::PHI];
    y[gr2::Weyl::UT] = -E/spt->get_metric()[gr2::Weyl::T][gr2::Weyl::T];
    y[gr2::Weyl::URHO] = 0.0;

    gr2::real norm2 = 0;
    for (int j = 0; j < 4; j++)
        norm2 += spt->get_metric()[j][j]*y[4+j]*y[4+j];
    ASSERT_GT(-1-norm2, 0);
    y[gr2::Weyl::UZ] = sqrtl((-1-norm2)/spt->get_metric()[gr2::Weyl::Z][gr2::Weyl::Z]);

    // initial conditions - steps
    const gr2::real dt = 0.7;
    int N = 10000;

    for (int i = 0; i < N; i++)
    {
        stepper->step(1, y, dt);
        spt->calculate_metric(y);
        // energy
        EXPECT_NEAR(E, -spt->get_metric()[gr2::Weyl::T][gr2::Weyl::T]*y[gr2::Weyl::UT], eps);
        // angular momentum
        EXPECT_NEAR(L, spt->get_metric()[gr2::Weyl::PHI][gr2::Weyl::PHI]*y[gr2::Weyl::UPHI], eps);
        // norm
        gr2::real norm2 = 0;
        for (int j = 0; j < 4; j++)
            norm2 += spt->get_metric()[j][j]*y[4+j]*y[4+j];
        EXPECT_NEAR(-1, norm2, eps);
    }
}

TEST(WeylSchwarzschildTrajectory, IntegralsOfMotionPlanarWithLambda)
{
    gr2::real eps = 1e-13;

    // prepare objects
    auto stepper = std::make_shared<gr2::DoPr853>();
    auto spt = std::make_shared<gr2::WeylSchwarzschild>(1.0, gr2::exact, gr2::diff);
    stepper->set_OdeSystem(spt);
    gr2::real y[9]{};

    // initial conditions - position
    y[gr2::Weyl::RHO] = sqrtl(16*14);
    
    // lambda
    spt->calculate_lambda_init(y);
    y[gr2::Weyl::LAMBDA] = spt->get_lambda();

    // initial conditions - velocity
    spt->calculate_metric(y);
    gr2::real L = 3.6823981191047921;
    y[gr2::Weyl::UPHI] = L/spt->get_metric()[gr2::Weyl::PHI][gr2::Weyl::PHI];
    y[gr2::Weyl::UT] = sqrtl((-1 - L*y[gr2::Weyl::UPHI])/spt->get_metric()[gr2::Weyl::T][gr2::Weyl::T]);
    gr2::real E = -spt->get_metric()[gr2::Weyl::T][gr2::Weyl::T]*y[gr2::Weyl::UT];
    ASSERT_NEAR(E, 0.9598686615055122147, 1e-15); // check if energy is correct

    // initial conditions - steps
    gr2::real dt = 0.2;
    int N = 5000;
    gr2::real r_min = 16;

    for (int i = 0; i < N; i++)
    {
        stepper->step(0, y, dt, false, nullptr, nullptr);
        spt->calculate_metric(y);
        // energy
        EXPECT_NEAR(E, -spt->get_metric()[gr2::Weyl::T][gr2::Weyl::T]*y[gr2::Weyl::UT], eps);
        // angular momentum
        EXPECT_NEAR(L, spt->get_metric()[gr2::Weyl::PHI][gr2::Weyl::PHI]*y[gr2::Weyl::UPHI], eps);// norm of four-velocity
        // norm2
        gr2::real norm2 = 0;
        for (int j = 0; j < 4; j++)
            norm2 += spt->get_metric()[j][j]*y[4+j]*y[4+j];
        EXPECT_NEAR(-1, norm2, eps);
        // lambda
        spt->calculate_lambda_init(y);
        EXPECT_NEAR(spt->get_lambda(), y[gr2::Weyl::LAMBDA], eps);

        r_min = std::min(r_min, y[gr2::Weyl::RHO]);
    }
    std::cout << r_min << std::endl;
}

TEST(WeylSchwarzschildTrajectory, IntegralsOfMotionGeneralWithLambda)
{
    gr2::real eps = 1e-13;

    // prepare objects
    auto stepper = std::make_shared<gr2::DoPr853>();
    auto spt = std::make_shared<gr2::WeylSchwarzschild>(1.0, gr2::exact, gr2::diff);
    stepper->set_OdeSystem(spt);
    gr2::real y[9]{};

    // initial conditions - position
    y[gr2::Weyl::RHO] = sqrtl(16*14);
    
    // lambda
    spt->calculate_lambda_init(y);
    y[gr2::Weyl::LAMBDA] = spt->get_lambda();

    // initial conditions - velocity
    spt->calculate_metric(y);
    gr2::real L = 3.6823981191047921;
    gr2::real E = 0.97;
    y[gr2::Weyl::UPHI] = L/spt->get_metric()[gr2::Weyl::PHI][gr2::Weyl::PHI];
    y[gr2::Weyl::UT] = -E/spt->get_metric()[gr2::Weyl::T][gr2::Weyl::T];
    y[gr2::Weyl::URHO] = 0.0;

    gr2::real norm2 = 0;
    for (int j = 0; j < 4; j++)
        norm2 += spt->get_metric()[j][j]*y[4+j]*y[4+j];
    ASSERT_GT(-1-norm2, 0);
    y[gr2::Weyl::UZ] = sqrtl((-1-norm2)/spt->get_metric()[gr2::Weyl::Z][gr2::Weyl::Z]);

    // initial conditions - steps
    gr2::real dt = 0.2;
    int N = 5000;

    for (int i = 0; i < N; i++)
    {
        stepper->step(0, y, dt, false, nullptr, nullptr);
        spt->calculate_metric(y);
        // energy
        EXPECT_NEAR(E, -spt->get_metric()[gr2::Weyl::T][gr2::Weyl::T]*y[gr2::Weyl::UT], eps);
        // angular momentum
        EXPECT_NEAR(L, spt->get_metric()[gr2::Weyl::PHI][gr2::Weyl::PHI]*y[gr2::Weyl::UPHI], eps);// norm of four-velocity
        // norm2
        gr2::real norm2 = 0;
        for (int j = 0; j < 4; j++)
            norm2 += spt->get_metric()[j][j]*y[4+j]*y[4+j];
        EXPECT_NEAR(-1, norm2, eps);
        // lambda
        spt->calculate_lambda_init(y);
        EXPECT_NEAR(spt->get_lambda(), y[gr2::Weyl::LAMBDA], eps);
    }
}

int main(int argc, char **argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}