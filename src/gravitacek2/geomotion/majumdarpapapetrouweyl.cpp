#include <cmath>
#include <stdexcept>
#include <iostream>

#include "gravitacek2/geomotion/majumadpapapetrouweyl.hpp"

namespace gr2
{
    MajumdarPapapetrouWeyl::MajumdarPapapetrouWeyl() : GeoMotion(4, 8)
    {

    }

    MajumdarPapapetrouWeyl::~MajumdarPapapetrouWeyl()
    {

    }

    real MajumdarPapapetrouWeyl::get_N_inv() const
    {
        return this->N_inv;
    }

    real MajumdarPapapetrouWeyl::get_N_inv_rho() const
    {
        return this->N_inv_rho;
    }

    real MajumdarPapapetrouWeyl::get_N_inv_z() const
    {
        return this->N_inv_z;
    }

    real MajumdarPapapetrouWeyl::get_N_inv_rhorho() const
    {
        return this->N_inv_rhorho;
    }

    real MajumdarPapapetrouWeyl::get_N_inv_rhoz() const
    {
        return this->N_inv_rhoz;
    }

    real MajumdarPapapetrouWeyl::get_N_inv_zz() const
    {
        return this->N_inv_zz;
    }

    void MajumdarPapapetrouWeyl::calculate_metric(const real *y)
    {
        if(!necessary_calculate(y, y_m, dim))
            return;

        real t = y[T];
        real phi = y[PHI];
        real rho = y[RHO];
        real z = y[Z];

        this->calculate_N_inv(y);
        real N = 1.0/N_inv;

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

        this->calculate_N_inv1(y);
        real N = 1.0/N_inv;

        christoffel_symbols[T][T][RHO] = -N_inv_rho*N;
        christoffel_symbols[T][RHO][T] = christoffel_symbols[T][T][RHO];
        christoffel_symbols[T][T][Z] = -N_inv_z*N;
        christoffel_symbols[T][Z][T] = christoffel_symbols[T][T][Z];
        christoffel_symbols[PHI][PHI][RHO] = N_inv_rho*N + 1.0/rho;
        christoffel_symbols[PHI][RHO][PHI] = christoffel_symbols[PHI][PHI][RHO];
        christoffel_symbols[PHI][PHI][Z] = N_inv_z*N;
        christoffel_symbols[PHI][Z][PHI] = christoffel_symbols[PHI][PHI][Z];
        christoffel_symbols[RHO][T][T] = -N_inv_rho*((N*N)*(N*N))*N;
        christoffel_symbols[RHO][PHI][PHI] = -rho*rho *N_inv_rho*N - rho;
        christoffel_symbols[RHO][RHO][RHO] = N_inv_rho*N;
        christoffel_symbols[RHO][RHO][Z] = N_inv_z*N;
        christoffel_symbols[RHO][Z][RHO] = christoffel_symbols[RHO][RHO][Z];
        christoffel_symbols[RHO][Z][Z] = -N_inv_rho*N;
        christoffel_symbols[Z][T][T] = -N_inv_z*((N*N)*(N*N))*N;
        christoffel_symbols[Z][PHI][PHI] = -rho*rho*N_inv_z*N;
        christoffel_symbols[Z][RHO][RHO] = -N_inv_z*N;
        christoffel_symbols[Z][RHO][Z] = N_inv_rho*N;
        christoffel_symbols[Z][Z][RHO] = christoffel_symbols[Z][RHO][Z];
        christoffel_symbols[Z][Z][Z] = N_inv_z*N;
    };

    void MajumdarPapapetrouWeyl::calculate_riemann_tensor(const real *y)
    {
        if(!necessary_calculate(y, y_r, dim))
            return;

        real t = y[T];
        real phi = y[PHI];
        real rho = y[RHO];
        real z = y[Z];

        this->calculate_N_inv2(y);
        real N = 1.0/N_inv;
        real N2 = N*N;
        real N6 = N*N*N*N*N*N;
        real rho_inv = 1.0/rho;

        riemann_tensor[T][PHI][T][PHI] = rho*(rho*N_inv_z*N_inv_z + (rho*N_inv_rho + N_inv)*N_inv_rho)*N2;
        riemann_tensor[T][PHI][PHI][T] = -riemann_tensor[T][PHI][T][PHI];
        riemann_tensor[T][RHO][T][RHO] = (N_inv*N_inv_rhorho - 3*N_inv_rho*N_inv_rho + N_inv_z*N_inv_z)*N2;
        riemann_tensor[T][RHO][RHO][T] = -riemann_tensor[T][RHO][T][RHO];
        riemann_tensor[T][RHO][T][Z] = (N_inv*N_inv_rhoz - 4*N_inv_rho*N_inv_z)*N2;
        riemann_tensor[T][RHO][Z][T] = -riemann_tensor[T][RHO][T][Z];
        riemann_tensor[T][Z][T][RHO] = riemann_tensor[T][RHO][T][Z];
        riemann_tensor[T][Z][RHO][T] = -riemann_tensor[T][Z][T][RHO];
        riemann_tensor[T][Z][T][Z] = (N_inv*N_inv_zz + N_inv_rho*N_inv_rho - 3*N_inv_z*N_inv_z)*N2;
        riemann_tensor[T][Z][Z][T] = -riemann_tensor[T][Z][T][Z];
        riemann_tensor[PHI][T][T][PHI] = (rho*N_inv_z*N_inv_z + (rho*N_inv_rho+N_inv)*N_inv_rho)*rho_inv*N6;
        riemann_tensor[PHI][T][PHI][T] = -riemann_tensor[PHI][T][T][PHI];
        riemann_tensor[PHI][RHO][PHI][RHO] = (-rho*N_inv*N_inv_rhorho+rho*N_inv_rho*N_inv_rho - rho*N_inv_z*N_inv_z - N_inv*N_inv_rho)*rho_inv*N2;
        riemann_tensor[PHI][RHO][RHO][PHI] = -riemann_tensor[PHI][RHO][PHI][RHO];
        riemann_tensor[PHI][RHO][PHI][Z] = (-N_inv*N_inv_rhoz + 2*N_inv_rho*N_inv_z)*N2;
        riemann_tensor[PHI][RHO][Z][PHI] = -riemann_tensor[PHI][RHO][PHI][Z]; 
        riemann_tensor[PHI][Z][PHI][RHO] = riemann_tensor[PHI][RHO][PHI][Z];
        riemann_tensor[PHI][Z][RHO][PHI] = -riemann_tensor[PHI][Z][PHI][RHO];
        riemann_tensor[PHI][Z][PHI][Z] = (-rho*N_inv*N_inv_zz + rho*N_inv_z*N_inv_z - (rho*N_inv_rho + N_inv)*N_inv_rho)*rho_inv*N2;
        riemann_tensor[PHI][Z][Z][PHI] = -riemann_tensor[PHI][Z][PHI][Z];
        riemann_tensor[RHO][T][T][RHO] = (N_inv*N_inv_rhorho - 3*N_inv_rho*N_inv_rho + N_inv_z*N_inv_z)*N6;
        riemann_tensor[RHO][T][RHO][T] = -riemann_tensor[RHO][T][T][RHO];
        riemann_tensor[RHO][T][T][Z] = (N_inv*N_inv_rhoz-4*N_inv_rho*N_inv_z)*N6;
        riemann_tensor[RHO][T][Z][T] = -riemann_tensor[RHO][T][T][Z];
        riemann_tensor[RHO][PHI][PHI][RHO] = (rho*N_inv*N_inv_rhorho - rho*N_inv_rho*N_inv_rho + rho*N_inv_z*N_inv_z + N_inv*N_inv_rho)*rho*N2;
        riemann_tensor[RHO][PHI][RHO][PHI] = -riemann_tensor[RHO][PHI][PHI][RHO];
        riemann_tensor[RHO][PHI][PHI][Z] = rho*rho*(N_inv*N_inv_rhoz - 2*N_inv_rho*N_inv_z)*N2;
        riemann_tensor[RHO][PHI][Z][PHI] = -riemann_tensor[RHO][PHI][PHI][Z];
        riemann_tensor[RHO][Z][RHO][Z] = (-(N_inv_rhorho + N_inv_zz)*N_inv+N_inv_rho*N_inv_rho+N_inv_z*N_inv_z)*N2;
        riemann_tensor[RHO][Z][Z][RHO] = -riemann_tensor[RHO][Z][RHO][Z];
        riemann_tensor[Z][T][T][RHO] = (N_inv*N_inv_rhoz-4*N_inv_rho*N_inv_z)*N6;
        riemann_tensor[Z][T][RHO][T] = -riemann_tensor[Z][T][T][RHO];
        riemann_tensor[Z][T][T][Z] = (N_inv*N_inv_zz + N_inv_rho*N_inv_rho - 3*N_inv_z*N_inv_z)*N6;
        riemann_tensor[Z][T][Z][T] = -riemann_tensor[Z][T][T][Z];
        riemann_tensor[Z][PHI][PHI][RHO] = (N_inv*N_inv_rhoz-2*N_inv_rho*N_inv_z)*rho*rho*N2;
        riemann_tensor[Z][PHI][RHO][PHI] = -riemann_tensor[Z][PHI][PHI][RHO];
        riemann_tensor[Z][PHI][PHI][Z] = (rho*N_inv*N_inv_zz-rho*N_inv_z*N_inv_z + (rho*N_inv_rho+N_inv)*N_inv_rho)*rho*N2;
        riemann_tensor[Z][PHI][Z][PHI] = -riemann_tensor[Z][PHI][PHI][Z];
        riemann_tensor[Z][RHO][RHO][Z] = ((N_inv_rhorho + N_inv_zz)*N_inv - N_inv_rho*N_inv_rho - N_inv_z*N_inv_z)*N2;
        riemann_tensor[Z][RHO][Z][RHO] = -riemann_tensor[Z][RHO][RHO][Z];
    }
}