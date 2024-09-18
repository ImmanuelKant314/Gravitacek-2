#include "gtest/gtest.h"
#include "gravitacek2/odesolver/ode.hpp"

TEST(odesolver_test, damped_harmonic_oscillator)
{
    class DampedHarmonicOscillator : public gr2::ODE
    {
        protected:
            gr2::REAL omega_0;
            gr2::REAL xi;

        public:
            DampedHarmonicOscillator(const gr2::REAL &omega_0, const gr2::REAL &xi):gr2::ODE(2) 
            {
                this->omega_0 = omega_0;
                this->xi = xi;
            }

            virtual void function(const gr2::REAL &t, const gr2::REAL y[], gr2::REAL dydt[]) override
            {
                dydt[0] = y[1];
                dydt[1] = 2*xi*omega_0*y[1] + omega_0*omega_0*y[0];
            }
    };

    DampedHarmonicOscillator osc = DampedHarmonicOscillator(1.5, 2);
    gr2::REAL y[] = {0.5, 1.5};
    gr2::REAL dydt[2];
    osc.function(0, y, dydt);

    gr2::REAL eps = 1e-15;
    ASSERT_LE(std::abs(dydt[0]-1.5), eps);
    ASSERT_LE(std::abs(dydt[1]-10.125), eps);
}

int main(int argc, char **argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}