#include "gravitacek2/setup.hpp"
#include "gravitacek2/odesolver/stepcontroller.hpp"

namespace gr2
{
    /**
     * @brief Step Controller based on relative and absolute error.
     * 
     * This controller is originally implemented in [GNU Scientific Library](https://www.gnu.org/software/gsl/doc/html/ode-initval.html#c.gsl_odeiv2_control_standard_new)
     * more details can be find there. Here about algoritm shortly. To calculate 
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
    class StandardStepController : public StepController
    {
    protected:
        int k;
        real eps_abs, eps_rel, a_y, a_dydt, S, factor;
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
        virtual real hadjust(const real y[], const real err[], const real dydt[], const real &h) override;
    };
} 
