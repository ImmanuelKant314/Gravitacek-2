#include <memory>
#include <stdexcept>

#include "gravitacek2/geomotion/weyl.hpp"

namespace gr2
{
    void CombinedWeyl::calculate_lambda_integral(const real *y)
    {
        calculate_lambda_from_inf_to_z(y, 1e-15);
    };

    CombinedWeyl::CombinedWeyl(std::vector<std::shared_ptr<Weyl>> sources):Weyl(gr2::integral, gr2::diff),sources(sources)
    {

    };

    CombinedWeyl::~CombinedWeyl()
    {

    };

    void CombinedWeyl::calculate_lambda_init(real const* y)
    {
        switch (this->lambda_eval_init)
        {
        case LambdaEvaluation::integral:
            this->calculate_lambda_integral(y);
            break;
        default:
            throw std::runtime_error("Calculating lambda this way is not possible");
            break;
        }
    };

    void CombinedWeyl::calculate_lambda_run(real const* y)
    {
        switch (this->lambda_eval_run)
        {
        case LambdaEvaluation::diff:
            this->calculate_lambda_diff(y);
            break;
        default:
            throw std::runtime_error("Calculating lambda this way is not possible");
            break;
        }
    };

    void CombinedWeyl::calculate_nu(const real* y)
    {
        this->nu = 0;
        for (auto s : this->sources)
        {
            s->calculate_nu(y);
            this->nu += s->get_nu();
        }
    };

    void CombinedWeyl::calculate_nu1(const real* y)
    {
        this->nu = 0;
        this->nu_rho = 0;
        this->nu_z = 0;
        for (auto s : this->sources)
        {
            s->calculate_nu1(y);
            this->nu += s->get_nu();
            this->nu_rho += s->get_nu_rho();
            this->nu_z += s->get_nu_z();
        }
    };


    void CombinedWeyl::calculate_nu2(const real* y)
    {
        this->nu = 0;
        this->nu_rho = 0;
        this->nu_z = 0;
        this->nu_rhorho = 0;
        this->nu_rhoz = 0;
        this->nu_zz = 0;
        for (auto s : this->sources)
        {
            s->calculate_nu2(y);
            this->nu += s->get_nu();
            this->nu_rho += s->get_nu_rho();
            this->nu_z += s->get_nu_z();
            this->nu_rhorho += s->get_nu_rhorho();
            this->nu_rhoz += s->get_nu_rhoz();
            this->nu_zz += s->get_nu_zz();
        }
    };
}