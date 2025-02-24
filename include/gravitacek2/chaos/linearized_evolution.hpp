#pragma once
#include "gravitacek2/geomotion/geomotion.hpp"

#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_cblas.h>

namespace gr2
{
    gsl_matrix* matrix_H(GeoMotion* spt, const real *y);

    gsl_matrix* time_corrected_matrix_H(GeoMotion* spt, const real *y);
}