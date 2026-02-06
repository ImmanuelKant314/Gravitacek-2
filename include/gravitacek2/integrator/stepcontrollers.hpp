/**
 * @file stepcontrollers.hpp
 * @author Karel Kraus
 * @brief Specific StepControllers for changing stepsize during integraion.
 * 
 * @copyright Copyright (c) 2026
 */

#pragma once
#include "gravitacek2/setup.hpp"
#include "gravitacek2/integrator/stepcontrollerbase.hpp"

namespace gr2
{
    /**
     * @brief Step Controller based on relative and absolute error.
     * 
     * This controller is originally implemented in [GNU Scientific
     * Library](https://www.gnu.org/software/gsl/doc/html/ode-initval.html#c.gsl_odeiv2_control_standard_new)
     * more details can be found there. Here about algoritm shortly. To
     * calculate new step size, we first calculate desired error for each ODE as
     * \f[
     * D_i = \epsilon_\text{abs} + \epsilon_\text{rel}(a_y \left|y_i\right| + a_{dydt}\left|y_i'\right|),
     * \f]
     * where \f$\epsilon_\text{abs}\f$, \f$\epsilon_\text{rel}\f$ are desired
     * relative and absolute error. Than we determine ratio \f$E_i/D_i\f$, where
     * \f$E_i = |y_\text{err}|\f$. We than find maximal value \f$E/D\f$.
     * 
     * If \f$E/D > 1.1\f$, than we calculate new step size as 
     * \f[
     * h_\text{new} = h_\text{old} S \left(\frac{E}{D}\right)^{-1/k},
     * \f]
     * 
     * where \f$S\f$ is a safety factor (we use 0.95) and \f$k\f$ is a
     * consintency order for a method for solving ODE. 
     * 
     * Else if \f$E/D < 0.5\f$, we calculate new step size as 
     * \f[
     * h_\text{new} = h_\text{old} S \left(\frac{E}{D}\right)^{-1/(k+1)}.
     * \f]
     * 
     * To avoid uncontrolable changes of step size, \f$h_\text{new}\f$ is
     * reduced to interval \f$[h_\text{old}/5, 5h_\text{old}]\f$.
     */
    class StandardStepController : public StepControllerBase
    {
    protected:
        int k;          //!<oder of stepper
        real eps_abs;   //!<desired absolute error
        real eps_rel;   //!<desired relative error
        real a_y;       //!<coefficient for `y`
        real a_dydt;    //!<coefficient for `dydt`
        real S;         //!<safety factor
        real factor;    //!<maximal relative change of step
        real *ratios;   //!<array for storing values \f$E_i/D_i\f$
    public:
        /**
         * @brief Construct a new Standard Step Controller object.
         * 
         * @param n number of ODE
         * @param k order of stepper
         * @param eps_abs desired absolute error
         * @param eps_rel desired relative error
         * @param a_y coefficient for `y`
         * @param a_dydt coefficient for `dydt`
         * @param S safety factor
         * @param factor maximal change factor
         */
        StandardStepController(const int &n, const int &k, const real &eps_abs, const real &eps_rel, const real &a_y, const real &a_dydt, const real &S = 0.95, const real &factor = 5);
        ~StandardStepController();
        virtual bool hadjust(const real y[], const real err[], const real dydt[], real &h) override;
    };

    /**
     * @brief StepController based on relative and absolute tolerance.
     * 
     * This controller is originally implemented in [Numerical
     * recipes](https://numerical.recipes/book.html), 
     * more details can be found there. We will here describe algorithm shortly. 
     * First of all we calculate scale for our error based on tolerance as
     * \f[
     *   \mbox{scale}_i = \mbox{atol} + \mbox{rtol} \left| y_i \right|.
     * \f]
     * Then we calculate error as
     * \f[
     *   \mbox{err} = \sqrt{\sum_{i}
     *   \left(\frac{\mbox{err}_i}{\mbox{scale}_i}\right)^2}.
     * \f]
     * If \f$\mbox{err}>1\f$, step is rejected, otherwise step is accepted. Next
     * step size can be calculated as
     * \f[
     * h_{\mbox{new}} = h_{\mbox{old}}S\left(\frac{1}{\mbox{err}}\right)^{1/k},
     * \f]
     * where \f$S\f$ is safety factor (because error estimate is not precise) 
     * and \f$k\f$ is order of steper used for calculation.
     */
    class StepControllerNR : public StepControllerBase
    {
    protected:
        int k;                  //!<order of stepper
        real atol;              //!<absolute error tolerance
        real rtol;              //!<relative error tolerance
        real S;                 //!<safety factor
        real factor_decrease;   //!<maximum decrease factor
        real factor_grow;       //!<maximum growth factor
        real err;               //!<variable for storing error
        real scale;             //!<variable for storing scale

    public:
        /**
         * @brief Construct a new StepControllerNR object  
         * 
         * @param n number of ordinary differential equations
         * @param k order of integrator
         * @param atol absolute tolerance of error
         * @param rtol relative tolerance of error
         * @param S safety factor, default 0.95
         * @param factor_decrease minimal growth factor, default 1/5
         * @param factor_grow maximal growth factor, default 10
         */
        StepControllerNR(const int &n, const int &k, const real &atol, const real &rtol, const real &S = 0.95, const real &factor_decrease = 1.0/5 ,const real &factor_grow = 10);
        // ~StepControllerNR();
        virtual bool hadjust(const real y[], const real err[], const real dydt[], real &h) override;
    };
} 

