#include <stdexcept>

#include "gravitacek2/chaos/linearized_evolution.hpp"

#include <gsl/gsl_linalg.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>

#include <stdio.h>
#include <gsl/gsl_matrix.h>

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

    real expected_growth(GeoMotion* spt, const real *y)
    {
        gsl_matrix *H;
        gsl_vector_complex *eig_vals;
        gsl_eigen_nonsymm_workspace *work_space;

        // calculate matrix H
        H = gr2::matrix_H(spt, y);
        gsl_matrix_transpose(H);
        eig_vals = gsl_vector_complex_alloc(8);
        work_space = gsl_eigen_nonsymm_alloc(8);

        // calculate eigen values
        int i = gsl_eigen_nonsymm(H, eig_vals, work_space);

        // find greates real part
        gr2::real value=0;
        for (int i = 0; i < work_space->n_evals; i++)
            value = std::max((double)value, GSL_REAL(gsl_vector_complex_get(eig_vals, i)));

        // free the matrix
        gsl_matrix_free(H);
        gsl_vector_complex_free(eig_vals);
        gsl_eigen_nonsymm_free(work_space);

        return value;
    }

    // gsl_matrix *time_corrected_matrix_H(GeoMotion *spt, const real *y)
    // {
    //     // constants
    //     int dim = spt->get_dim();
    //     int n = 2*dim;

    //     // arrays and special constants
    //     gsl_matrix *H_sym, *H, *fw; // matrix H
    //     gsl_matrix *f, *w; // vectors
    //     gsl_matrix *A, *B, *g, *gH, *gH_sym; // help matrix
    //     real f_norm2 = 0;
    //     real fHf = 0;
        
    //     try
    //     {
    //         // calculate H
    //         H = matrix_H(spt, y);

    //         // TODO: recalculate f (this should be (
    //         // calculate f (done)
    //         spt->calculate_christoffel_symbols(y);
    //         f = gsl_matrix_alloc(n, 1);
    //         for (int i = 0; i< dim; i++)
    //         {
    //             gsl_matrix_set(f, i, 0, y[dim+i]);
    //             gr2::real value = 0;
    //             for (int j = 0; j < dim; j++)
    //                 for (int k = 0; k < dim; k++)
    //                     value += -spt->get_christoffel_symbols()[i][j][k]*y[dim+j]*y[dim+k];
    //             gsl_matrix_set(f, dim+i, 0, value);
    //         }

    //         // calculate g (done)
    //         spt->calculate_metric(y);
    //         g = gsl_matrix_alloc(n, n);
    //         gsl_matrix_set_zero(g);
    //         for (int i = 0; i < dim; i++)
    //             for (int j = 0; j < dim; j++)
    //             {
    //                 gsl_matrix_set(g, i, j, spt->get_metric()[i][j]);
    //                 gsl_matrix_set(g, dim+i, dim+j, spt->get_metric()[i][j]);
    //             }

    //         // calculate f_norm2 (done)
    //         for (int i = 0; i < n; i++)
    //             for (int j = 0; j < n; j++)
    //                 f_norm2 += gsl_matrix_get(g, i, j)*gsl_matrix_get(f, i, 0)*gsl_matrix_get(f, j, 0);

    //         // calculate gH (done)
    //         gH = gsl_matrix_alloc(n, n);
    //         gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, g, H, 0.0, gH);
            
    //         // calculate gH symmetrize (done)
    //         gH_sym = gsl_matrix_alloc(n, n);
    //         gsl_matrix_transpose_memcpy(gH_sym, gH);
    //         gsl_matrix_add(gH_sym, gH);

    //         // calculate A (done)
    //         A = gsl_matrix_alloc(n, 1);
    //         gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, -1.0/(f_norm2), gH_sym, f, 0.0, A);

    //         // calculate fgHf (done)
    //         for (int i = 0; i < n; i++)
    //             for (int j = 0; j < n; j++)
    //                 fHf += gsl_matrix_get(gH, i, j)*gsl_matrix_get(f, i, 0)*gsl_matrix_get(f, j, 0);

    //         // calculate B (done)
    //         B = gsl_matrix_alloc(n, 1);
    //         gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, fHf/(f_norm2*f_norm2), g, f, 0.0, B);

    //         // calculate w
    //         w = gsl_matrix_alloc(n, 1);
    //         gsl_matrix_memcpy(w, A);
    //         gsl_matrix_add(w, B);

    //         // calculate fw
    //         fw =  gsl_matrix_alloc(n, n);
    //         gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, f, w, 0, fw);

    //         // calculate H
    //         gsl_matrix_add(H, fw);

    //         // free memory
    //         gsl_matrix_free(f);
    //         gsl_matrix_free(g);
    //         gsl_matrix_free(gH);
    //         gsl_matrix_free(gH_sym);
    //         gsl_matrix_free(A);
    //         gsl_matrix_free(B);
    //         gsl_matrix_free(w);
    //         gsl_matrix_free(fw);

    //         // return
    //         return H;
    //     }
    //     catch(const std::exception& e)
    //     {
    //         // free memory
    //         gsl_matrix_free(H);
    //         gsl_matrix_free(f);
    //         gsl_matrix_free(g);
    //         gsl_matrix_free(gH);
    //         gsl_matrix_free(gH_sym);
    //         gsl_matrix_free(A);
    //         gsl_matrix_free(B);
    //         gsl_matrix_free(w);
    //         gsl_matrix_free(fw);
    //         throw e;
    //     }   
    // };

    gr2::real max_norm_growth(GeoMotion *spt, const real *y)
    {
        // constants
        int dim = spt->get_dim();
        int n = 2*dim;

        // arrays and special constants
        gsl_matrix *H, *gH, *g, *gHT;   // matrix H
        gsl_matrix *mod_g;              // modified metric
        gsl_matrix *P;                  // projector
        gsl_matrix *matrix;             // final matrix
        gsl_matrix *tmp;                // temporary matrix for multiplication
        gr2::real u_up[4]{}, u_down[4]{};

        gsl_vector *eig_vals;
        gsl_eigen_symm_workspace *work_space;

        try
        {
            // calculate H
            H = matrix_H(spt, y);

            tmp = gsl_matrix_alloc(n, n);

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

            // sum gH and (gH).T to gH
            gsl_matrix_add(gH, gHT);

            // TODO: Calculate velocities
            for (int i = 0; i < dim; i++)
                u_up[i] = y[i+dim];

            for (int i = 0; i < dim; i++)
            {
                u_down[i] = 0;
                for (int j = 0; j < dim; j++)
                    u_down[i] += spt->get_metric()[i][j]*u_up[j];
            }

            // TODO: Calculate projector
            P = gsl_matrix_alloc(n, n);
            gsl_matrix_set_zero(P);
            for (int i = 0; i< dim; i++)
                for (int j = 0; j< dim; j++)
                {
                    if (i == j)
                    {
                        gsl_matrix_set(P, i, j, 1+u_up[i]*u_down[j]);
                        gsl_matrix_set(P, i+dim, j+dim, 1+u_up[i]*u_down[j]);
                    }
                    else
                    {
                        gsl_matrix_set(P, i, j, u_up[i]*u_down[j]);
                        gsl_matrix_set(P, i+dim, j+dim, u_up[i]*u_down[j]);
                    }
                }

            // TODO: Calculate modified metric
            mod_g = gsl_matrix_alloc(n, n);
            gsl_matrix_set_zero(mod_g);
            for (int i = 0; i< dim; i++)
                for (int j = 0; j< dim; j++)
                {
                    gsl_matrix_set(mod_g, i, j, spt->get_metric()[i][j] + 2*u_down[i]*u_down[j]);
                    gsl_matrix_set(mod_g, i+4, j+4, spt->get_metric()[i][j] + 2*u_down[i]*u_down[j]);
                }

            // TODO: Calculate Cholesky decomposition
            gsl_linalg_cholesky_decomp1(mod_g);

            // delete upper diagonal
            for (int i = 0; i < n; i++)
                for (int j = i+1; j < n; j++)
                    gsl_matrix_set(mod_g, i, j, 0);

            // TODO: Calculate inverse of L
            gsl_linalg_tri_invert(CblasLower, CblasNonUnit, mod_g);

            // save P*L^-1 to P
            gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, P, mod_g, 0.0, tmp);
            gsl_matrix_memcpy(P, tmp);

            // TODO:Calculate resulting matrix
            matrix = gsl_matrix_alloc(n, n);
            gsl_matrix_set_identity(matrix);
            gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, P, gH, 0.0, tmp);
            gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, tmp, P, 0.0, matrix);

            // TODO: calculate value
            // calculate matrix H
            eig_vals = gsl_vector_alloc(8);
            work_space = gsl_eigen_symm_alloc(8);

            // calculate eigen values
            int status = gsl_eigen_symm(matrix, eig_vals, work_space);
            gr2::real value=0;
            for (int i = 0; i < 8; i++)
                value = std::max((double)value, gsl_vector_get(eig_vals, i));

            // delete matrices
            gsl_matrix_free(H);
            gsl_matrix_free(gH);
            gsl_matrix_free(g);
            gsl_matrix_free(gHT);
            gsl_matrix_free(mod_g);
            gsl_matrix_free(P);
            gsl_matrix_free(tmp);
            gsl_matrix_free(matrix);
            gsl_vector_free(eig_vals);
            gsl_eigen_symm_free(work_space);

            // return
            return value*0.5;
        }
        catch(const std::exception& e)
        {
            gsl_matrix_free(H);
            gsl_matrix_free(gH);
            gsl_matrix_free(g);
            gsl_matrix_free(gHT);
            gsl_matrix_free(mod_g);
            gsl_matrix_free(P);
            gsl_matrix_free(tmp);
            gsl_matrix_free(matrix);
            gsl_vector_free(eig_vals);
            gsl_eigen_symm_free(work_space);
            throw e;
        }
        
    }
}