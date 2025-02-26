#include <cmath>

#include "gtest/gtest.h"
#include "gravitacek2/integrator/odesystem.hpp"
#include "gravitacek2/integrator/event.hpp"
#include "gravitacek2/integrator/stepperbase.hpp"
#include "gravitacek2/integrator/steppers.hpp"

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

void linear_regression(gr2::real x_data[], gr2::real y_data[], int N, gr2::real &a, gr2::real &b)
{
    gr2::real sum_xy = 0, sum_x = 0, sum_y = 0, sum_x2 = 0;
    for (int i = 0; i < N; i++)
    {
        sum_x += x_data[i];
        sum_y += y_data[i];
        sum_xy += x_data[i]*y_data[i];
        sum_x2 += x_data[i]*x_data[i];
    }
    a = (N*sum_xy - sum_x*sum_y)/(N*sum_x2 - sum_x*sum_x);
    b = (sum_y*sum_x2 - sum_x*sum_xy)/(N*sum_x2-sum_x*sum_x);
}

class StepperTestCase
{
public:
    std::string name;                           // name of the integrator
    std::shared_ptr<gr2::StepperBase> stepper;  // tested stepper
    gr2::real min_exp, max_exp;                 // range of exponents 10**x for steps size for measuring OrderOfIntegrator
    gr2::real eps_integration;                  // maximal error during integration
    gr2::real eps_ord;                          // maximal error for orders
    gr2::real h;                                // stepsize for integration

    StepperTestCase(std::string name, std::shared_ptr<gr2::StepperBase> stepper, gr2::real min_exp, gr2::real max_exp, gr2::real eps_integration, gr2::real eps_ord, gr2::real h):
    name(name), stepper(stepper), min_exp(min_exp), max_exp(max_exp), eps_integration(eps_integration), eps_ord(eps_ord), h(h)
    {}
};

class GeneralStepperTest:public testing::TestWithParam<StepperTestCase>
{

};

TEST_P(GeneralStepperTest, OrderOfIntegrator)
{
    // parameters and variables
    gr2::real omega0 = 1.5, xi = 1.0;
    gr2::real x0 = 0.5, v0 = 1.5;
    gr2::real dydt[2];

    int N = 10;
    gr2::real min_exp = GetParam().min_exp, max_exp = GetParam().max_exp;
    gr2::real x_data[N], y_data[N];
    gr2::real exp, h, err;
    gr2::real total_sum_orders = 0;
    gr2::real y[2];

    // ODE
    auto osc = std::make_shared<DampedHarmonicOscillator>(omega0, xi);

    // Stepper
    auto stepper = GetParam().stepper;
    stepper->set_OdeSystem(osc);

    // calculate errors
    for (int i = 0; i < N; i++)
    {
        exp = min_exp + (max_exp-min_exp)/(N-1)*i;
        h = powl(10, exp);
        y[0] = x0;
        y[1] = v0;

        stepper->step(0, y, h);
        err = fabs(y[0] - exactDampedHarmonicOscillator(h, omega0, xi, x0, v0));
        x_data[i] = exp;
        y_data[i] = log10(err);
    }
    gr2::real order, bias;
    linear_regression(x_data, y_data, N, order, bias);
    EXPECT_NEAR(order, stepper->get_order()+1, GetParam().eps_ord);
}

TEST_P(GeneralStepperTest, OrderOfError)
{
    // parameters and variables
    gr2::real omega0 = 1.5, xi = 1.0;
    gr2::real x0 = 0.5, v0 = 1.5;
    gr2::real dydt[2];

    int N = 10;
    gr2::real min_exp = GetParam().min_exp, max_exp = GetParam().max_exp;
    gr2::real x_data[N], y_data[N];
    gr2::real exp, h, err;
    gr2::real total_sum_orders = 0;
    gr2::real y[2], error[2];

    // ODE
    auto osc = std::make_shared<DampedHarmonicOscillator>(omega0, xi);

    // Stepper
    auto stepper = GetParam().stepper;
    stepper->set_OdeSystem(osc);

    // calculate errors
    for (int i = 0; i < N; i++)
    {
        exp = min_exp + (max_exp-min_exp)/(N-1)*i;
        h = powl(10, exp);
        y[0] = x0;
        y[1] = v0;

        stepper->step_err(0, y, h, error);
        x_data[i] = exp;
        y_data[i] = log10(std::fabs(error[0]));
    }
    gr2::real order, bias;
    linear_regression(x_data, y_data, N, order, bias);
    EXPECT_NEAR(order, stepper->get_err_order(), GetParam().eps_ord);
}

// correctness of the integration
TEST_P(GeneralStepperTest, Integration)
{
    // parameters and variables
    gr2::real omega0 = 1.5, xi = 1.0;
    gr2::real x0 = 0.5, v0 = 1.5;
    gr2::real h = GetParam().h;
    gr2::real y[2]{x0, v0}, dydt[2];

    // ODE
    auto osc = std::make_shared<DampedHarmonicOscillator>(omega0, xi);

    // Stepper
    auto stepper = GetParam().stepper;
    stepper->set_OdeSystem(osc);

    // test steps
    for (int i = 0; i < int(7/h); i++)
    {
        ASSERT_NEAR(y[0], exactDampedHarmonicOscillator(h*i, omega0, xi, x0, v0), GetParam().eps_integration);
        stepper->step(i*h, y, h);
    }
}

// correctness of the integration (with error)
TEST_P(GeneralStepperTest, IntegrationWithError)
{
    // parameters and variables
    gr2::real omega0 = 1.5, xi = 1.0;
    gr2::real x0 = 0.5, v0 = 1.5;
    gr2::real h = GetParam().h;
    gr2::real y[2]{x0, v0}, dydt[2], err[2];

    // ODE
    auto osc = std::make_shared<DampedHarmonicOscillator>(omega0, xi);

    // Stepper
    auto stepper = GetParam().stepper;
    stepper->set_OdeSystem(osc);

    // test steps
    for (int i = 0; i < int(7/h); i++)
    {
        ASSERT_NEAR(y[0], exactDampedHarmonicOscillator(h*i, omega0, xi, x0, v0), GetParam().eps_integration);
        stepper->step_err(i*h, y, h, err);
    }
}

// test of dense output
// TEST_P(GeneralStepperTest, DenseOutput)
// {

// }

void PrintTo(const StepperTestCase& testcase, std::ostream* os) {
    *os << "0";
}

std::string MyTestNameGenerator(const ::testing::TestParamInfo<StepperTestCase>& info) {
    return info.param.name;
};

auto test_cases = testing::Values(
    StepperTestCase("RK4", std::make_shared<gr2::RK4>(false), -3, -1, 1e-7, 0.25, 0.002),
    StepperTestCase("DoPr853", std::make_shared<gr2::DoPr853>(), -1, 1, 1e-7, 0.6, 0.01)
);

INSTANTIATE_TEST_SUITE_P(
    StepperTest,
    GeneralStepperTest,
    test_cases,
    MyTestNameGenerator
);

int main(int argc, char **argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}