#include <stdexcept>

#include "gravitacek2/chaos/linearized_evolution.hpp"

namespace gr2
{
    gsl_matrix *matrix_H(GeoMotion* spt, const real *y)
    {
        gsl_matrix *matrix;
        try
        {
            int dim = spt->get_dim();
            int n = 2*dim;

            matrix = gsl_matrix_alloc(n, n);
            spt->calculate_riemann_tensor(y);

            // make array zeros
            gsl_matrix_set_zero(matrix);
            
            // add ones
            for (int i = 0; i < dim; i++)
                gsl_matrix_set(matrix, i, i+4, 1);

            // add riemann tensor
            for (int i = 0; i < dim; i++)
                for (int j = 0; j < dim; j++)
                {
                    gr2::real value = 0;
                    for (int k = 0; k < dim; k++)
                        for (int l = 0; l < dim; l++)
                            value += -spt->get_riemann_tensor()[i][k][j][l]*y[dim+k]*y[dim+l];
                    gsl_matrix_set(matrix, dim+i, j, value);
            }
        }
        catch(const std::exception& e)
        {
            gsl_matrix_free(matrix);
            throw e;
        }

        return matrix;
    };

    gsl_matrix *time_corrected_matrix_H(GeoMotion *spt, const real *y)
    {
        // constants
        int dim = spt->get_dim();
        int n = 2*dim;

        // arrays and special constants
        gsl_matrix *H_sym, *H, *fw; // matrix H
        gsl_matrix *f, *w; // vectors
        gsl_matrix *A, *B, *g, *gH, *gH_sym; // help matrix
        real f_norm2 = 0;
        real fHf = 0;
        
        try
        {
            // calculate H
            H = matrix_H(spt, y);

            // calculate f (done)
            spt->calculate_christoffel_symbols(y);
            f = gsl_matrix_alloc(n, 1);
            for (int i = 0; i< dim; i++)
            {
                gsl_matrix_set(f, i, 0, y[dim+i]);
                gr2::real value = 0;
                for (int j = 0; j < dim; j++)
                    for (int k = 0; k < dim; k++)
                        value += -spt->get_christoffel_symbols()[i][j][k]*y[dim+j]*y[dim+k];
                gsl_matrix_set(f, dim+i, 0, value);
            }

            // calculate g (done)
            spt->calculate_metric(y);
            g = gsl_matrix_alloc(n, n);
            gsl_matrix_set_zero(g);
            for (int i = 0; i < dim; i++)
                for (int j = 0; j < dim; j++)
                {
                    gsl_matrix_set(g, i, j, spt->get_metric()[i][j]);
                    gsl_matrix_set(g, dim+i, dim+j, spt->get_metric()[i][j]);
                }

            // calculate f_norm2 (done)
            for (int i = 0; i < n; i++)
                for (int j = 0; j < n; j++)
                    f_norm2 += gsl_matrix_get(g, i, j)*gsl_matrix_get(f, i, 0)*gsl_matrix_get(f, j, 0);

            // calculate gH (done)
            gH = gsl_matrix_alloc(n, n);
            gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, g, H, 0.0, gH);
            
            // calculate gH symmetrize (done)
            gH_sym = gsl_matrix_alloc(n, n);
            gsl_matrix_transpose_memcpy(gH_sym, gH);
            gsl_matrix_add(gH_sym, gH);

            // calculate A (done)
            A = gsl_matrix_alloc(n, 1);
            gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, -1.0/(f_norm2), gH_sym, f, 0.0, A);

            // calculate fgHf (done)
            for (int i = 0; i < n; i++)
                for (int j = 0; j < n; j++)
                    fHf += gsl_matrix_get(gH, i, j)*gsl_matrix_get(f, i, 0)*gsl_matrix_get(f, j, 0);

            // calculate B (done)
            B = gsl_matrix_alloc(n, 1);
            gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, fHf/(f_norm2*f_norm2), g, f, 0.0, B);

            // calculate w
            w = gsl_matrix_alloc(n, 1);
            gsl_matrix_memcpy(w, A);
            gsl_matrix_add(w, B);

            // calculate fw
            fw =  gsl_matrix_alloc(n, n);
            gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, f, w, 0, fw);

            // calculate H
            gsl_matrix_add(H, fw);

            // free memory
            gsl_matrix_free(f);
            gsl_matrix_free(g);
            gsl_matrix_free(gH);
            gsl_matrix_free(gH_sym);
            gsl_matrix_free(A);
            gsl_matrix_free(B);
            gsl_matrix_free(w);
            gsl_matrix_free(fw);

            // return
            return H;
        }
        catch(const std::exception& e)
        {
            // free memory
            gsl_matrix_free(H);
            gsl_matrix_free(f);
            gsl_matrix_free(g);
            gsl_matrix_free(gH);
            gsl_matrix_free(gH_sym);
            gsl_matrix_free(A);
            gsl_matrix_free(B);
            gsl_matrix_free(w);
            gsl_matrix_free(fw);
            throw e;
        }   
    };

    gsl_matrix *norm_growth(GeoMotion *spt, const real *y)
    {
        // constants
        int dim = spt->get_dim();
        int n = 2*dim;

        // arrays and special constants
        gsl_matrix *H, *gH, *g, *gHT; // matrix H

        try
        {
            // calculate H
            H = time_corrected_matrix_H(spt, y);

            // calculate g (done)
            spt->calculate_metric(y);
            g = gsl_matrix_alloc(n, n);
            gsl_matrix_set_zero(g);
            for (int i = 0; i < dim; i++)
                for (int j = 0; j < dim; j++)
                {
                    gsl_matrix_set(g, i, j, spt->get_metric()[i][j]);
                    gsl_matrix_set(g, dim+i, dim+j, spt->get_metric()[i][j]);
                }

            // calculate gH
            gH = gsl_matrix_alloc(n, n);
            gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, g, H, 0, gH);

            // calculate (gH).T
            gHT = gsl_matrix_alloc(n, n);
            gsl_matrix_transpose_memcpy(gHT, gH);
        
            // calculate final result
            gsl_matrix_add(gH, gHT);
            gsl_matrix_scale(gH, 0.5);

            // delete matrices
            gsl_matrix_free(H);
            gsl_matrix_free(g);
            gsl_matrix_free(gHT);

            // return
            return gH;
        }
        catch(const std::exception& e)
        {
            gsl_matrix_free(H);
            gsl_matrix_free(gH);
            gsl_matrix_free(g);
            gsl_matrix_free(gHT);
            throw e;
        }
        
    }
}