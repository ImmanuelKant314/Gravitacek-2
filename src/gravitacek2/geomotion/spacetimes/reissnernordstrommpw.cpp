#include <stdexcept>
#include <cmath>

#include "gravitacek2/mymath.hpp"
#include "gravitacek2/geomotion/spacetimes.hpp"

namespace gr2
{
    ReissnerNordstromMPW::ReissnerNordstromMPW(const real& M)
    {
        this->M = M;
    };

    ReissnerNordstromMPW::~ReissnerNordstromMPW()
    {

    };

    void ReissnerNordstromMPW::calculate_N_inv(const real* y) 
    {
        // coordinates
        real rho = y[RHO];
        real z = y[Z];
        real d = sqrtl(rho*rho + z*z);
        real d_inv = 1.0/d;

        // lapse
        this->N_inv = 1.0 + M*d_inv;
    };

    void ReissnerNordstromMPW::calculate_N_inv1(const real* y)
    {
        // coordinates
        real rho = y[RHO];
        real z = y[Z];
        real d = sqrtl(rho*rho + z*z);
        real d_inv = 1.0/d;
        real M_d_inv = 1.0/(M+d);

        // lapse
        this->N_inv = 1.0 + M*d_inv;

        // derivative of lapse
        this->N_inv_rho = -M*rho*d_inv*d_inv*d_inv;
        this->N_inv_z = -M*z*d_inv*d_inv*d_inv;
    };

    void ReissnerNordstromMPW::calculate_N_inv2(const real* y)
    {
        // coordinates
        real rho = y[RHO];
        real z = y[Z];
        real d = sqrtl(rho*rho + z*z);
        real d_inv = 1.0/d;
        real M_d_inv = 1.0/(M+d);
        real d_inv5 = d_inv*d_inv*d_inv*d_inv*d_inv;

        // lapse
        this->N_inv = 1.0 + M*d_inv;

        // derivative of lapse
        this->N_inv_rho = -M*rho*d_inv*d_inv*d_inv;
        this->N_inv_z = -M*z*d_inv*d_inv*d_inv;

        // second derivative of lapse
        this->N_inv_rhorho = M*(2*rho*rho-z*z)*d_inv5;
        this->N_inv_rhoz = 3*M*rho*z*d_inv5;
        this->N_inv_zz = M*(2*z*z-rho*rho)*d_inv5;
    };
}