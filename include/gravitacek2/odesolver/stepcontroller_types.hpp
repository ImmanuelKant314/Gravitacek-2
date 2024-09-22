#include "gravitacek2/setup.hpp"
#include "gravitacek2/odesolver/stepcontroller.hpp"

namespace gr2
{
    class StandardStepController : public StepController
    {
    protected:
        int k;
        real eps_abs, eps_rel, a_y, a_dydt, S, factor;
        real *ratios;
    public:
        StandardStepController(const int &n, const int &k, const real &eps_abs, const real &eps_rel, const real &a_y, const real &a_dydt, const real &S = 0.95, const real &factor = 5);
        ~StandardStepController();
        virtual real hadjust(const real y[], const real err[], const real dydt[], const real &h) override;
    };
} 
