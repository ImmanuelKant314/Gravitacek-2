#pragma once
#include "gravitacek2/odesolver/stepper.hpp"

namespace gr2
{
    /**
     * @brief Stepper using algorithm RK4.
     * 
     * <table>
     *   <caption>Buther's tableau for RK4</caption>
     *   <tr>
     *     <td>0</td> <td>0</td> <td>0</td> <td>0</td> <td>0</td> 
     *   </tr>
     *   <tr>
     *     <td>1/2</td> <td>1/2</td> <td>0</td> <td>0</td> <td>0</td>
     *   </tr>
     *   <tr>
     *     <td>1/2</td> <td>0</td> <td>1/2</td> <td>0</td> <td>0</td> 
     *   </tr>
     *   <tr>
     *     <td>1</td> <td>0</td> <td>0</td> <td>1</td> <td>0</td> 
     *   </tr>
     *   <tr>
     *     <td></td> <td>1/6</td> <td>1/3</td> <td>1/3</td> <td>1/6</td>
     *   </tr>
     * </table>
     * 
     */
    class RK4 : public Stepper
    {
    protected:
        REAL *k1, *k2, *k3, *k4;
    public:
        RK4();
        ~RK4();
        virtual void set_ODE(ODE& ode);
        virtual void reset();
        virtual void step(const REAL &t, REAL y[], const REAL &h, const REAL dydt_in[] = nullptr, REAL dydt_out[] = nullptr) override;
        virtual int get_order() const override;
    };
}