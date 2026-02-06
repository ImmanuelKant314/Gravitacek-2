/**
 * @file odesystems.hpp
 * @author Karel Kraus
 * @brief Classes representing specific ODEs.
 * 
 * @copyright Copyright (c) 2026
 * 
 */

#pragma once
#include "gravitacek2/integrator/odesystem.hpp"

namespace gr2
{
    /**
     * @brief OdeSystem class for damped harmonic oscilator.
     * 
     * Damped harmonic oscilator is represented by differential equation
     * \f[
     * \dv{t} \begin{pmatrix} x \\ v \end{pmatrix} = 
     * \begin{pmatrix}
     * v \\ -2\xi v - \omega_0^2x
     * \end{pmatrix},
     * \f]
     * where \f$\omega_0\f$ is undamped angular frequency of the oscillator and
     * \f$\xi\f$ is damping ratio.
     */
    class DampedHarmonicOscillator : public gr2::OdeSystem
    {
        protected:
            gr2::real omega0;   //!<undamped angular frequency
            gr2::real xi;       //!<damping ratio

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