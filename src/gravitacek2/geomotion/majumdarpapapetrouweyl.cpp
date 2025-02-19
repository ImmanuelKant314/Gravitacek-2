#include <cmath>
#include <stdexcept>

#include "gravitacek2/geomotion/majumadpapapetrouweyl.hpp"

namespace gr2
{
    MajumdarPapapetrouWeyl::MajumdarPapapetrouWeyl() : GeoMotion(4, 8)
    {

    }

    MajumdarPapapetrouWeyl::~MajumdarPapapetrouWeyl()
    {

    }

    real MajumdarPapapetrouWeyl::get_N() const
    {
        return this->N;
    }

    real MajumdarPapapetrouWeyl::get_invN() const
    {
        return this->N_inv;
    }

    real MajumdarPapapetrouWeyl::get_N_rho() const
    {
        return this->N_rho;
    }

    real MajumdarPapapetrouWeyl::get_N_z() const
    {
        return this->N_z;
    }

    real MajumdarPapapetrouWeyl::get_N_rhorho() const
    {
        return this->N_rhorho;
    }

    real MajumdarPapapetrouWeyl::get_N_rhoz() const
    {
        return this->N_rhoz;
    }

    real MajumdarPapapetrouWeyl::get_N_zz() const
    {
        return this->N_zz;
    }

    void MajumdarPapapetrouWeyl::calculate_metric(const real *y)
    {
        if(!necessary_calculate(y, y_m, dim))
            return;

        real t = y[T];
        real phi = y[PHI];
        real rho = y[RHO];
        real z = y[Z];

        this->calculate_N(y);
        
        metric[T][T] = -N*N;
        metric[PHI][PHI] = (rho*rho)*(N_inv*N_inv);
        metric[RHO][RHO] = N_inv*N_inv;
        metric[Z][Z] = metric[RHO][RHO];
    }

    void MajumdarPapapetrouWeyl::calculate_christoffel_symbols(const real *y)
    {
        if(!necessary_calculate(y, y_c, dim))
            return;
        
        real t = y[T];
        real phi = y[PHI];
        real rho = y[RHO];
        real z = y[Z];

        this->calculate_N1(y);

        christoffel_symbols[T][T][RHO] = N_rho*N_inv;
        christoffel_symbols[T][RHO][T] = christoffel_symbols[T][T][RHO];
        christoffel_symbols[T][T][Z] = N_z*N_inv;
        christoffel_symbols[T][Z][T] = christoffel_symbols[T][T][Z];
        christoffel_symbols[PHI][PHI][RHO] = -N_rho*N_inv + 1.0/rho;
        christoffel_symbols[PHI][RHO][PHI] = christoffel_symbols[PHI][PHI][RHO];
        christoffel_symbols[PHI][PHI][Z] = -N_z*N_inv;
        christoffel_symbols[PHI][Z][PHI] = christoffel_symbols[PHI][PHI][Z];
        christoffel_symbols[RHO][T][T] = N*N*N*N_rho;
        christoffel_symbols[RHO][PHI][PHI] = rho*rho *N_rho*N_inv - rho;
        christoffel_symbols[RHO][RHO][RHO] = -N_rho*N_inv;
        christoffel_symbols[RHO][RHO][Z] = -N_z*N_inv;
        christoffel_symbols[RHO][Z][RHO] = christoffel_symbols[RHO][RHO][Z];
        christoffel_symbols[RHO][Z][Z] = N_rho*N_inv;
        christoffel_symbols[Z][T][T] = N*N*N*N_z;
        christoffel_symbols[Z][PHI][PHI] = rho*rho*N_z*N_inv;
        christoffel_symbols[Z][RHO][RHO] = N_z*N_inv;
        christoffel_symbols[Z][RHO][Z] = -N_rho*N_inv;
        christoffel_symbols[Z][Z][RHO] = christoffel_symbols[Z][RHO][Z];
        christoffel_symbols[Z][Z][Z] = -N_z*N_inv;
    };

    void MajumdarPapapetrouWeyl::calculate_riemann_tensor(const real *y)
    {
        if(!necessary_calculate(y, y_r, dim))
            return;

        real t = y[T];
        real phi = y[PHI];
        real rho = y[RHO];
        real z = y[Z];

        this->calculate_N2(y);
        real N_inv2 = N_inv*N_inv;
        real rho_inv = 1.0/rho;

        riemann_tensor[T][PHI][T][PHI] = rho*(rho*N_z*N_z + (rho*N_rho - N)*N_rho)*N_inv2;
        riemann_tensor[T][PHI][PHI][T] = -riemann_tensor[T][PHI][T][PHI];
        riemann_tensor[T][RHO][T][RHO] = (-N*N_rhorho - N_rho*N_rho + N_z*N_z)*N_inv2;
        riemann_tensor[T][RHO][RHO][T] = -riemann_tensor[T][RHO][T][RHO];
        riemann_tensor[T][RHO][T][Z] = -(N*N_rhoz + 2*N_rho*N_z)*N_inv2;
        riemann_tensor[T][RHO][Z][T] = -riemann_tensor[T][RHO][T][Z];
        riemann_tensor[T][Z][T][RHO] = riemann_tensor[T][RHO][T][Z];
        riemann_tensor[T][Z][RHO][T] = -riemann_tensor[T][Z][T][RHO];
        riemann_tensor[T][Z][T][Z] = (-N*N_zz + N_rho*N_rho - N_z*N_z)*N_inv2;
        riemann_tensor[T][Z][Z][T] = -riemann_tensor[T][Z][T][Z];
        riemann_tensor[PHI][T][T][PHI] = (rho*N_z*N_z + (rho*N_rho-N)*N_rho)*N*N*rho_inv;
        riemann_tensor[PHI][T][PHI][T] = -riemann_tensor[PHI][T][T][PHI];
        riemann_tensor[PHI][RHO][PHI][RHO] = (rho*N*N_rhorho-rho*N_rho*N_rho - rho*N_z*N_z + N*N_rho)*rho_inv*N_inv2;
        riemann_tensor[PHI][RHO][RHO][PHI] = -riemann_tensor[PHI][RHO][PHI][RHO];
        riemann_tensor[PHI][RHO][PHI][Z] = N_rhoz*N_inv;
        riemann_tensor[PHI][RHO][Z][PHI] = -riemann_tensor[PHI][RHO][PHI][Z]; 
        riemann_tensor[PHI][Z][PHI][RHO] = riemann_tensor[PHI][RHO][PHI][Z];
        riemann_tensor[PHI][Z][RHO][PHI] = -riemann_tensor[PHI][Z][PHI][RHO];
        riemann_tensor[PHI][Z][PHI][Z] = (rho*N*N_zz - rho*N_z*N_z - (rho*N_rho - N)*N_rho)*rho_inv*N_inv2;
        riemann_tensor[PHI][Z][Z][PHI] = -riemann_tensor[PHI][Z][PHI][Z];
        riemann_tensor[RHO][T][T][RHO] = (-N*N_rhorho - N_rho*N_rho + N_z*N_z)*N*N;
        riemann_tensor[RHO][T][RHO][T] = -riemann_tensor[RHO][T][T][RHO];
        riemann_tensor[RHO][T][T][Z] = (-N*N_rhoz-2*N_rho*N_z)*N*N;
        riemann_tensor[RHO][T][Z][T] = -riemann_tensor[RHO][T][T][Z];
        riemann_tensor[RHO][PHI][PHI][RHO] = (-rho*N*N_rhorho + rho*N_rho*N_rho + rho*N_z*N_z - N*N_rho)*rho*N_inv2;
        riemann_tensor[RHO][PHI][RHO][PHI] = -riemann_tensor[RHO][PHI][PHI][RHO];
        riemann_tensor[RHO][PHI][PHI][Z] = -(rho*rho*N_rhoz)*N_inv;
        riemann_tensor[RHO][PHI][Z][PHI] = -riemann_tensor[RHO][PHI][PHI][Z];
        riemann_tensor[RHO][Z][RHO][Z] = ((N_rhorho + N_zz)*N-N_rho*N_rho-N_z*N_z)*N_inv2;
        riemann_tensor[RHO][Z][Z][RHO] = -riemann_tensor[RHO][Z][RHO][Z];
        riemann_tensor[Z][T][T][RHO] = (-N*N_rhoz-2*N_rho*N_z)*N*N;
        riemann_tensor[Z][T][RHO][T] = -riemann_tensor[Z][T][T][RHO];
        riemann_tensor[Z][T][T][Z] = (-N*N_zz + N_rho*N_rho - N_z*N_z)*N*N;
        riemann_tensor[Z][T][Z][T] = -riemann_tensor[Z][T][T][Z];
        riemann_tensor[Z][PHI][PHI][RHO] = -rho*rho*N_rhoz*N_inv;
        riemann_tensor[Z][PHI][RHO][PHI] = -riemann_tensor[Z][PHI][PHI][RHO];
        riemann_tensor[Z][PHI][PHI][Z] = (-rho*N*N_zz+rho*N_z*N_z + (rho*N_rho-N)*N_rho)*rho*N_inv2;
        riemann_tensor[Z][PHI][Z][PHI] = -riemann_tensor[Z][PHI][PHI][Z];
        riemann_tensor[Z][RHO][RHO][Z] = (-(N_rhorho + N_zz)*N + N_rho*N_rho + N_z*N_z)*N_inv2;
        riemann_tensor[Z][RHO][Z][RHO] = -riemann_tensor[Z][RHO][RHO][Z];
    }
}