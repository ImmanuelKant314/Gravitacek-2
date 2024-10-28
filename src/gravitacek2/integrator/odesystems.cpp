#include "gravitacek2/integrator/odesystems.hpp"

DampedHarmonicOscillator::DampedHarmonicOscillator(const gr2::real &omega0, const gr2::real &xi):gr2::OdeSystem(2) 
{
    this->omega0 = omega0;
    this->xi = xi;
};

void DampedHarmonicOscillator::function(const gr2::real &t, const gr2::real y[], gr2::real dydt[])
{
    dydt[0] = y[1];
    dydt[1] = -2*this->xi*y[1] - this->omega0*this->omega0*y[0];
};