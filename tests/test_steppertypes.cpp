#include <cmath>

#include "gtest/gtest.h"
#include "gravitacek2/odesolver/ode.hpp"
#include "gravitacek2/odesolver/event.hpp"
#include "gravitacek2/odesolver/stepper.hpp"
#include "gravitacek2/odesolver/steper_types.hpp"

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

TEST(OrderOfIntegrator, RK4)
{
    // parameters and variables
    gr2::real omega0 = 1.5, xi = 1.0;
    gr2::real x0 = 0.5, v0 = 1.5;
    gr2::real dydt[2];
    gr2::real eps = 0.25;

    gr2::real min_exp = -3, max_exp = -1;
    int N = 10;
    gr2::real x_data[N], y_data[N];
    gr2::real exp, h, err;
    gr2::real total_sum_orders = 0;
    gr2::real y[2];

    // ODE
    DampedHarmonicOscillator osc = DampedHarmonicOscillator(omega0, xi);

    // Stepper
    gr2::RK4 stepper = gr2::RK4();
    stepper.set_ODE(osc);

    // calculate errors
    for (int i = 0; i < N; i++)
    {
        exp = min_exp + (max_exp-min_exp)/(N-1)*i;
        h = powl(10, exp);
        y[0] = x0;
        y[1] = v0;

        stepper.step(0, y, h);
        err = fabs(y[0] - exactDampedHarmonicOscillator(h, omega0, xi, x0, v0));
        x_data[i] = exp;
        y_data[i] = log10(err);
    }
    gr2::real order, bias;
    linear_regression(x_data, y_data, N, order, bias);
    EXPECT_NEAR(order, stepper.get_order()+1, eps);
}

TEST(OrderOfIntegrator, RK4Error)
{
    // parameters and variables
    gr2::real omega0 = 1.5, xi = 1.0;
    gr2::real x0 = 0.5, v0 = 1.5;
    gr2::real dydt[2];
    gr2::real eps = 0.25;

    gr2::real min_exp = -3, max_exp = -1;
    int N = 10;
    gr2::real x_data[N], y_data[N];
    gr2::real exp, h, err;
    gr2::real total_sum_orders = 0;
    gr2::real y[2], error[2];

    // ODE
    DampedHarmonicOscillator osc = DampedHarmonicOscillator(omega0, xi);

    // Stepper
    gr2::RK4 stepper = gr2::RK4();
    stepper.set_ODE(osc);

    // calculate errors
    for (int i = 0; i < N; i++)
    {
        exp = min_exp + (max_exp-min_exp)/(N-1)*i;
        h = powl(10, exp);
        y[0] = x0;
        y[1] = v0;

        stepper.step_err(0, y, h, error);
        err = fabs(y[0] - exactDampedHarmonicOscillator(h, omega0, xi, x0, v0));
        x_data[i] = exp;
        y_data[i] = log10(err);
    }
    gr2::real order, bias;
    linear_regression(x_data, y_data, N, order, bias);
    EXPECT_NEAR(order, stepper.get_order()+1, eps);
}

int main(int argc, char **argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}