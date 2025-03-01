#pragma once
#include "gravitacek2/integrator/odesystem.hpp"

namespace gr2
{
    class DampedHarmonicOscillator : public gr2::OdeSystem
    {
        protected:
            gr2::real omega0;
            gr2::real xi;

        public:
            /**
             * @brief Construct a new Damped Harmonic Oscillator.
             * 
             * @param omega0 undamped angular frequency of the oscillator
             * @param xi damping ratio
             */
            DampedHarmonicOscillator(const gr2::real &omega0, const gr2::real &xi);

            virtual void function(const gr2::real &t, const gr2::real y[], gr2::real dydt[]); 
    };
}