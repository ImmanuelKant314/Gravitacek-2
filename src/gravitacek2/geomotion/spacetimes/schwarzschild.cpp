#include "gravitacek2/geomotion/spacetimes.hpp"

#include <cmath>

namespace gr2
{
    Schwarzschild::Schwarzschild(real M) : GeoMotion(4, 8)
    {
        this->M = M;
    }

    void Schwarzschild::calculate_metric(const real *y)
    {
        real t = y[T];
        real r = y[R];
        real theta = y[THETA];
        real phi = y[PHI];

        real sin_val = sinl(theta);

        metric[T][T] = -(1 - 2 * M / r);
        metric[R][R] = 1.0 / (1 - 2.0 * M / r);
        metric[THETA][THETA] = r * r;
        metric[PHI][PHI] = metric[THETA][THETA]*sin_val*sin_val;
    }

    void Schwarzschild::calculate_christoffel_symbols(const real *y)
    {
        real t = y[T];
        real r = y[R];
        real theta = y[THETA];
        real phi = y[PHI];

        real sin_val = sinl(theta);
        real cos_val = cosl(theta);

        // ========== Christoffels symbols =========
        christoffel_symbols[T][T][R] = M / (r*(r - 2 * M));
        christoffel_symbols[T][R][T] = christoffel_symbols[T][T][R];
        christoffel_symbols[R][T][T] = M*(r-2*M)/(r*r*r);
        christoffel_symbols[R][R][R] = M/(r*(2*M-r));
        christoffel_symbols[R][THETA][THETA] = 2*M-r;
        christoffel_symbols[R][PHI][PHI] = (2*M-r)*sin_val*sin_val;
        christoffel_symbols[THETA][R][THETA] = 1.0/r;
        christoffel_symbols[THETA][THETA][R] = christoffel_symbols[THETA][R][THETA];
        christoffel_symbols[THETA][PHI][PHI] = - sin_val*cos_val;
        christoffel_symbols[PHI][R][PHI] = christoffel_symbols[THETA][R][THETA];
        christoffel_symbols[PHI][PHI][R] = christoffel_symbols[PHI][R][PHI];
        christoffel_symbols[PHI][THETA][PHI] = cos_val/sin_val;
        christoffel_symbols[PHI][PHI][THETA] = christoffel_symbols[PHI][THETA][PHI];
    }

    void Schwarzschild::calculate_riemann_tensor(const real *y)
    {
        real t = y[T];
        real r = y[R];
        real theta = y[THETA];
        real phi = y[PHI];

        real sin_val = sinl(theta);
        real cos_val = cosl(theta);

        // ========== Riemann tensor ========== 
        riemann_tensor[T][R][T][R] = 2*M/(r*r*(r-2*M));
        riemann_tensor[T][R][R][T] = -riemann_tensor[T][R][T][R];
        riemann_tensor[T][THETA][T][THETA] = -M/r;
        riemann_tensor[T][THETA][THETA][T] = -riemann_tensor[T][THETA][T][THETA];
        riemann_tensor[T][PHI][T][PHI] = - M*sin_val*sin_val/r;
        riemann_tensor[T][PHI][PHI][T] = -riemann_tensor[T][PHI][T][PHI];
        riemann_tensor[R][T][T][R] = 2*M*(r-2*M)/(r*r*r*r);
        riemann_tensor[R][T][R][T] = -riemann_tensor[R][T][T][R];
        riemann_tensor[R][THETA][R][THETA] = - M/r;
        riemann_tensor[R][THETA][THETA][R] = -riemann_tensor[R][THETA][R][THETA];
        riemann_tensor[R][PHI][R][PHI] = - M*sin_val*sin_val/r;
        riemann_tensor[R][PHI][PHI][R] = -riemann_tensor[R][PHI][R][PHI];
        riemann_tensor[THETA][T][T][THETA] = M*(2*M-r)/(r*r*r*r);
        riemann_tensor[THETA][T][THETA][T] = -riemann_tensor[THETA][T][T][THETA];
        riemann_tensor[THETA][R][R][THETA] = M/(r*r*(r-2*M));
        riemann_tensor[THETA][R][THETA][R] = -riemann_tensor[THETA][R][R][THETA];
        riemann_tensor[THETA][PHI][THETA][PHI] = 2*M*sin_val*sin_val/r;
        riemann_tensor[THETA][PHI][PHI][THETA] = -riemann_tensor[THETA][PHI][THETA][PHI];
        riemann_tensor[PHI][T][T][PHI] = M*(2*M-r)/(r*r*r*r);
        riemann_tensor[PHI][T][PHI][T] = -riemann_tensor[PHI][T][T][PHI];
        riemann_tensor[PHI][R][R][PHI] = M/(r*r*(r-2*M));
        riemann_tensor[PHI][R][PHI][R] = -riemann_tensor[PHI][R][R][PHI];
        riemann_tensor[PHI][THETA][THETA][PHI] = - 2*M/r;
        riemann_tensor[PHI][THETA][PHI][THETA] = -riemann_tensor[PHI][THETA][THETA][PHI];
    }
}