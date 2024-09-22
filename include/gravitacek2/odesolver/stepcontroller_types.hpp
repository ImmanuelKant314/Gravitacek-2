#include "gravitacek2/setup.hpp"
#include "gravitacek2/odesolver/stepcontroller.hpp"

namespace gr2
{
    class StandardStepController : public StepController
    {
    protected:
        int k;
        REAL eps_abs, eps_rel, a_y, a_dydt, S, factor;
        REAL *ratios;
    public:
        StandardStepController(const int &n, const int &k, const REAL &eps_abs, const REAL &eps_rel, const REAL &a_y, const REAL &a_dydt, const REAL &S = 0.95, const REAL &factor = 5);
        ~StandardStepController();
        virtual REAL hadjust(const REAL y[], const REAL err[], const REAL dydt[], const REAL &h) override;
    };
} 
