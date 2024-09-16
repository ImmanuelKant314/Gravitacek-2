#include "ode.hpp"

namespace gr2
{
    /**
     * @brief Construct a new ODE object
     *
     * Default constructor takes number of equations. 
     */
    ODE::ODE(const int &n)
    {
        this->n = n;
    }

    /**
     * @brief Return number of equations.
     *
     * @return int number of equations
     */
    int ODE::get_n() const
    {
        return this->n;
    }
}