#pragma once
#include "gravitacek2/setup.hpp"
#include "gravitacek2/integrator/stepcontrollerbase.hpp"

namespace gr2
{
    /**
     * @brief Step Controller based on relative and absolute error.
     * 
     * This controller is originally implemented in [GNU Scientific Library](https://www.gnu.org/software/gsl/doc/html/ode-initval.html#c.gsl_odeiv2_control_standard_new)
     * more details can be found there. Here about algoritm shortly. To calculate 
     * new step size, we first calculate desired error for each ODE as
     * \f[
     * D_i = \epsilon_{abs} + \epsilon_{rel}(a_y \left|y_i\right| + a_{dydt}\left|y_i'\right|)
     * \f]
     * Than we determine ratio \f$E_i/D_i\f$ and find maximal value \f$E/D\f$.
     * 
     * If \f$E/D > 1.1\f$, than we calculate new step size as 
     * \f[
     *   h_{new} = h_{old} S \left(\frac{E}{D}\right)^{-1/k},
     * \f]

     * where \f$S\f$ is a safety factor of 0.95 and \f$k\f$ is a consintency 
     * order for a method for solving ODE. 
     * 
     * Else if \f$E/D < 0.5\f$, we calculate new step size as 
     * \f[
     *   h_{new} = h_{old} S \left(\frac{E}{D}\right)^{-1/(k+1)}.
     * \f]
     * 
     * To avoid uncontrolable changes of step size, \f$h_{new}\f$ is reduced
     * to interval \f$[h_{old}/5, 5h_{old}]\f$
     */
    class StandardStepController : public StepControllerBase
    {
    protected:
        int k;
        real eps_abs;
        real eps_rel;
        real a_y;
        real a_dydt;
        real S;
        real factor;
        real *ratios;
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
     * @brief Step Controller based on relative and absolute tolerance.
     * 
     * This controller is originally implemented in [Numerical recipes](https://numerical.recipes/book.html), 
     * more details can be found there. We will here describe algorithm shortly. 
     * First of all we calculate scale for our error based on tolerance as
     * \f[
     *   \mbox{scale}_i = \mbox{atol} + \mbox{rtol} \left| y_i \right|.
     * \f]
     * Then we calculate error as
     * \f[
     *   \mbox{err} = \sqrt{\sum_{i} \left(\frac{\mbox{err}_i}{\mbox{scale}_i}\right)^2}.
     * \f]
     * If \f$\mbox{err}>1\f$, step is rejected, otherwise step is accepted. Next step 
     * size can be calculated as
     * \f[
     * h_{\mbox{new}} = h_{\mbox{old}}S\left(\frac{1}{\mbox{err}}\right)^{1/k},
     * \f]
     * where \f$S\f$ is safety factor (because error estimate is not precise) 
     * and \f$k\f$ is order of steper used for calculation.
     */
    class StepControllerNR : public StepControllerBase
    {
    protected:
        real atol;              //!<absolute error tolerance
        real rtol;              //!<relative error tolerance
        int k;                  //!<order of stepper
        real S;                 //!<safety factor
        real factor_decrease;   //!<maximum decrease factor
        real factor_grow;       //!<maximum growth factor
        real err;               //!<variable for calculating error
        real scale;             //!<variable for calculating scale

    public:
        /**
         * @brief Construct a new Step Controller NR object  
         * 
         * @param n number of ordinary differential equations
         * @param k order of error
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

