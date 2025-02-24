#include "gtest/gtest.h"
#include "gravitacek2/setup.hpp"
#include "gravitacek2/chaos/linearized_evolution.hpp"
#include "gravitacek2/geomotion/spacetimes.hpp"

#include <gsl/gsl_linalg.h>

TEST(MatrixOfLinearizedEvolution, WeylSchwarzschild)
{
    gr2::real eps = 1e-10;
    gr2::real value;

    gr2::WeylSchwarzschild spt(1.0);
    gr2::real y[9] = {}; 
    y[gr2::Weyl::RHO] = 10;
    y[gr2::Weyl::Z] = 0;
    y[gr2::Weyl::UT] = 1.190472574611857e+00;
    y[gr2::Weyl::UPHI] = 3.071259328415933e-02;
    y[gr2::Weyl::URHO] = 1.000000000000000e-01;
    y[gr2::Weyl::UZ] = 1.663404206336647e-01;

    gr2::real H_test[4][4] = {
        {9.223591795963500e-5,-3.308860773972330e-3,2.133391827840170e-4,-1.774346470096790e-4},
        {2.219468244087360e-5,-9.009281293604130e-4,-2.751932170531810e-6,9.155151096031570e-6},
        {-1.445314775963320e-4,-2.779451492237140e-4,1.830764066209640e-3,-1.490455561880790e-5},
        {1.202071338908940e-4,9.246702606991880e-4,-1.490455561880790e-5,-1.022071854808860e-3}
    };

    // initialize lambda
    spt.calculate_lambda_init(y);
    y[gr2::Weyl::LAMBDA] = spt.get_lambda();

    // calculate matrix
    gsl_matrix* matrix = gr2::matrix_H(&spt, y);

    // check equality
    for (int i = 0; i < 4; i++)
        for (int j = 0; j < 4; j++)
        {
           EXPECT_NEAR(gsl_matrix_get(matrix, i, j), 0, eps);
           EXPECT_NEAR(gsl_matrix_get(matrix, i+4, j+4), 0, eps);
           value = (int)(i == j);
           EXPECT_NEAR(gsl_matrix_get(matrix, i, j+4), value, eps + std::labs(value)*eps);
           value = H_test[i][j];
           EXPECT_NEAR(gsl_matrix_get(matrix, i+4, j), value, eps + std::labs(value)*eps);
        }

    // delete
    gsl_matrix_free(matrix);
}

TEST(MatrixOfLinearizedEvolution, MathematicalCondition)
{
    gr2::real eps = 1e-10;
    gr2::real value;

    gr2::WeylSchwarzschild spt(1.0);
    gr2::real y[9] = {}; 
    y[gr2::Weyl::RHO] = 10;
    y[gr2::Weyl::Z] = 0;
    y[gr2::Weyl::UT] = 1.190472574611857e+00;
    y[gr2::Weyl::UPHI] = 3.071259328415933e-02;
    y[gr2::Weyl::URHO] = 1.000000000000000e-01;
    y[gr2::Weyl::UZ] = 1.663404206336647e-01;
    
    // initialize lambda
    spt.calculate_lambda_init(y);
    y[gr2::Weyl::LAMBDA] = spt.get_lambda();

    // calculate matrix
    gsl_matrix* matrix = gr2::time_corrected_matrix_H(&spt, y);
    
    // calculate acceleration
    gr2::real dydt[8]{};
    spt.calculate_christoffel_symbols(y);
    for (int i = 0; i < 4; i++)
    {
        dydt[i] = y[i+4];
        for (int j = 0; j < 4; j++)
            for (int k = 0; k < 4; k++)
                dydt[i+4] += -spt.get_christoffel_symbols()[i][j][k]*y[j+4]*y[k+4];
    }

    // calculate total metric
    gr2::real G[8][8]{};
    for (int i = 0; i<4; i++)
        for (int j = 0; j<4; j++)
        {
            G[i][j] = spt.get_metric()[i][j];
            G[i+4][j+4] = spt.get_metric()[i][j];
        }

    // calculate gH
    gr2::real gH[8][8]{};
    spt.calculate_metric(y);
    for (int i = 0; i<8; i++)
        for (int j = 0; j<8; j++)
            for (int k = 0; k<8; k++)
                gH[i][j] += G[i][k]*gsl_matrix_get(matrix, k, j);

    // check condition
    for (int i = 0; i<8; i++)
    {
        value=0;
        for (int j = 0; j<8; j++)
            value += (gH[i][j] + gH[j][i])*dydt[j];
        EXPECT_NEAR(value, 0, eps);
    }

    // delete
    gsl_matrix_free(matrix);
}

int main(int argc, char **argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}