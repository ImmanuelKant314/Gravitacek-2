#pragma once

/**
 * @brief Functions and object for Gravitacek 2.
 * 
*/
namespace gr2
{
// ========== Types ========== 

typedef long double real;               //!<default type for real numbers
typedef real (*realfunction)(real);     //!<type of real function


// ========== Constants ========== 

const real pi   = 3.1415926535897932385;    //!<value of \f$\pi\f$
const real pi_2 = 1.5707963267948966192;    //!<value of \f$\pi/2\f$
const real pi_4 = 0.78539816339744830962;   //!<value of \f$\pi/4\f$
const real e    = 2.7182818284590452354;    //!<value of \f$e\f$
}
