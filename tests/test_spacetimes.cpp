#include <iostream>
#include <fstream>

#include "gtest/gtest.h"
#include "gravitacek2/setup.hpp"
#include "gravitacek2/geomotion/spacetimes.hpp"

// ==================== Schwarzschild ==================== 

gr2::real y_schwarzschild[4] = {};
gr2::real metric_schwarzschild[4][4] = {};
gr2::real christoffel_symbols_schwarzschild[4][4][4] = {};
gr2::real riemann_tensor_schwarzschild[4][4][4][4] = {};
gr2::real eps = 1e-15;
std::string line;

TEST(Schwarzschild, Metric)
{
    gr2::Schwarzschild spacetime = gr2::Schwarzschild(0.3);
    spacetime.calculate_metric(y_schwarzschild);
    for (int i = 0; i < 4; i++)
        for (int j = 0; j < 4; j++)
            EXPECT_NEAR(spacetime.get_metric()[i][j], metric_schwarzschild[i][j], eps);
}

TEST(Schwarzschild, ChristoffelSymbols)
{
    gr2::Schwarzschild spacetime = gr2::Schwarzschild(0.3);
    spacetime.calculate_christoffel_symbols(y_schwarzschild);
    for (int i = 0; i < 4; i++)
        for (int j = 0; j < 4; j++)
            for(int k = 0; k < 4; k++)
                EXPECT_NEAR(spacetime.get_christoffel_symbols()[i][j][k], christoffel_symbols_schwarzschild[i][j][k], eps);
}

TEST(Schwarzschild, RiemannTensor)
{
    gr2::Schwarzschild spacetime = gr2::Schwarzschild(0.3);
    spacetime.calculate_riemann_tensor(y_schwarzschild);
    for (int i = 0; i < 4; i++)
        for (int j = 0; j < 4; j++)
            for (int k = 0; k < 4; k++)
                for (int l = 0; l < 4; l++)
                    EXPECT_NEAR(spacetime.get_riemann_tensor()[i][j][k][l], riemann_tensor_schwarzschild[i][j][k][l], eps);
}

// ==================== WeylSchwarzschild ==================== 

gr2::real y_weylschwarzschild[4] = {};
gr2::real metric_weylschwarzschild[4][4] = {};
gr2::real christoffel_symbols_weylschwarzschild[4][4][4] = {};
gr2::real riemann_tensor_weylschwarzschild[4][4][4][4] = {};

TEST(WeylSchwarzschild, Potencial)
{
    gr2::WeylSchwarzschild spacetime = gr2::WeylSchwarzschild(0.3);
    spacetime.calculate_nu(y_weylschwarzschild);
    EXPECT_NEAR(spacetime.get_nu(), -0.029413667983693322, eps);
}

TEST(WeylSchwarzchild, DerivativesOfPotential)
{
    gr2::WeylSchwarzschild spacetime = gr2::WeylSchwarzschild(0.3);
    spacetime.calculate_nu1(y_weylschwarzschild);
    EXPECT_NEAR(spacetime.get_nu(), -0.029413667983693322, eps);
    EXPECT_NEAR(spacetime.get_nu_rho(), 0.0028276099493778893, eps);
    EXPECT_NEAR(spacetime.get_nu_z(), 0.0005650330027458911, eps);
}

TEST(WeylSchwarzschild, LambdaExact)
{
    gr2::WeylSchwarzschild spacetime = gr2::WeylSchwarzschild(0.3);
    spacetime.calculate_lambda_init(y_weylschwarzschild);
    EXPECT_NEAR(spacetime.get_lambda(), -4.159049449726778e-04, eps);
}

TEST(WeylSchwarzschild, LambdaIntegral)
{
    gr2::WeylSchwarzschild spacetime = gr2::WeylSchwarzschild(0.3, gr2::integral);
    spacetime.calculate_lambda_init(y_weylschwarzschild);
    EXPECT_NEAR(spacetime.get_lambda(), -0.000415904944973, eps);
}

TEST(WeylSchwarzschild, Metric)
{
    gr2::WeylSchwarzschild spacetime = gr2::WeylSchwarzschild(0.3,gr2::exact,gr2::exact);
    spacetime.calculate_metric(y_weylschwarzschild);
    for (int i = 0; i < 4; i++)
        for (int j = 0; j < 4; j++)
            EXPECT_NEAR(spacetime.get_metric()[i][j], metric_weylschwarzschild[i][j], eps);
}

TEST(WeylSchwarzschild, ChristoffelSymbols)
{
    gr2::WeylSchwarzschild spacetime = gr2::WeylSchwarzschild(0.3,gr2::exact, gr2::exact);
    spacetime.calculate_christoffel_symbols(y_weylschwarzschild);
    for (int i = 0; i < 4; i++)
        for (int j = 0; j < 4; j++)
            for(int k = 0; k < 4; k++)
                EXPECT_NEAR(spacetime.get_christoffel_symbols()[i][j][k], christoffel_symbols_weylschwarzschild[i][j][k], 4*eps);
}

//TODO: WeylSchwarzschild - Riemannův tensor

// ==================== BachWeylRing ==================== 
gr2::real y_bachweylring[4] = {};

TEST(BachWeylRing, Potential)
{
    gr2::real eps = 1e-15;
    gr2::BachWeylRing spacetime = gr2::BachWeylRing(0.3, 20);
    spacetime.calculate_nu(y_bachweylring);
    EXPECT_NEAR(spacetime.get_nu(), -0.015000654984792, eps);
}

TEST(BachWeylRing, DerivativesOfPotential)
{
    gr2::real eps = 1e-15;
    gr2::BachWeylRing spacetime = gr2::BachWeylRing(0.3, 20);
    spacetime.calculate_nu1(y_bachweylring);
    EXPECT_NEAR(spacetime.get_nu(), -0.015000654984792, eps);
    EXPECT_NEAR(spacetime.get_nu_rho(), -0.000009372089805, eps);
    EXPECT_NEAR(spacetime.get_nu_z(), 0.000011262025116, eps);
}

// ==================== inverted Kuzmin-Toomre disks ==================== 
TEST(InvertedKuzminToomreDisk, Potential)
{
    gr2::real eps = 1e-15;
    gr2::real M = 3.000000000000000e-01;
    gr2::real b = 5.000000000000000e+00;

    gr2::real y[4] = {};
    y[2] = 1.000000000000000e+01;
    y[3] = 7.000000000000000e+00;

    gr2::real nu[] = {
        -1.684419608953263e-02,
        -1.526271395131956e-02,
        -1.407264603121147e-02};

    for (int n = 1; n <=3; n++)
    {
        gr2::InvertedKuzminToomreDisk spacetime = gr2::InvertedKuzminToomreDisk(n, M, b);
        spacetime.calculate_nu(y);
        EXPECT_NEAR(spacetime.get_nu(), nu[n-1], eps);
    }
}

TEST(InvertedKuzminToomreDisk, DerivativesOfPotential)
{
    gr2::real eps = 1e-15;
    gr2::real M = 3.000000000000000e-01;
    gr2::real b = 5.000000000000000e+00;

    gr2::real y[4] = {};
    y[2] = 1.000000000000000e+01;
    y[3] = 7.000000000000000e+00;

    gr2::real nu[] = {
        -1.684419608953263e-02,
        -1.526271395131956e-02,
        -1.407264603121147e-02};

    gr2::real nu_rho[] = {
        4.967837823249628e-04,
        3.380353453900308e-04,
        2.380667922163150e-04};

    gr2::real nu_z[] = {
        7.929185304901022e-04,
        6.774218538243885e-04,
        5.863221910056083e-04};

    for (int n = 1; n <=3; n++)
    {
        gr2::InvertedKuzminToomreDisk spacetime = gr2::InvertedKuzminToomreDisk(n, M, b);
        spacetime.calculate_nu1(y);
        EXPECT_NEAR(spacetime.get_nu(), nu[n-1], eps);
        EXPECT_NEAR(spacetime.get_nu_rho(), nu_rho[n-1], eps);
        EXPECT_NEAR(spacetime.get_nu_z(), nu_z[n-1], eps);
    }
}

TEST(InvertedKuzminToomreDisk, Lambda)
{

}

// ==================== inverted Morgan-Morgan disks ==================== 
TEST(InvertedMorganMorganDisk, Potential)
{
    gr2::real eps = 1e-15;
    gr2::real M = 3.000000000000000e-01;
    gr2::real b = 5.000000000000000e+00;

    gr2::real y[4] = {};
    y[2] = 1.000000000000000e+01;
    y[3] = 7.000000000000000e+00;

    gr2::real nu[] = {
        -1.747474699849155e-02,
        -1.557819025295148e-02,
        -1.423530400664234e-02};

    for (int n = 1; n <=3; n++)
    {
        gr2::InvertedMorganMorganDisk spacetime = gr2::InvertedMorganMorganDisk(n, M, b);
        spacetime.calculate_nu(y);
        EXPECT_NEAR(spacetime.get_nu(), nu[n-1], eps);
    }
}

TEST(InvertedMorganMorganDisk, DerivativesOfPotential)
{
    gr2::real eps = 1e-15;
    gr2::real M = 3.000000000000000e-01;
    gr2::real b = 5.000000000000000e+00;

    gr2::real y[4] = {};
    y[2] = 1.000000000000000e+01;
    y[3] = 7.000000000000000e+00;

    gr2::real nu[] = {
        -1.747474699849155e-02,
        -1.557819025295148e-02,
        -1.423530400664234e-02};

    gr2::real nu_rho[] = {
        4.816239271321822e-04,
        3.006630640524549e-04,
        1.977496695848732e-04};

    gr2::real nu_z[] = {
        8.677415301569396e-04,
        7.121903757523738e-04,
        6.000699761341162e-04};

    for (int n = 1; n <=3; n++)
    {
        gr2::InvertedMorganMorganDisk spacetime = gr2::InvertedMorganMorganDisk(n, M, b);
        spacetime.calculate_nu1(y);
        EXPECT_NEAR(spacetime.get_nu(), nu[n-1], eps);
        EXPECT_NEAR(spacetime.get_nu_rho(), nu_rho[n-1], eps);
        EXPECT_NEAR(spacetime.get_nu_z(), nu_z[n-1], eps);
    }
}

int main(int argc, char **argv)
{
    // ========== Prepare schwarzschild ========== 
    // coordinates
    y_schwarzschild[0] = 1.000000000000000;
    y_schwarzschild[1] = 4.000000000000000;
    y_schwarzschild[2] = 0.500000000000000;
    y_schwarzschild[3] = 0.300000000000000;

    // metric
    metric_schwarzschild[0][0] = -0.850000000000000;
    metric_schwarzschild[1][1] = 1.176470588235294;
    metric_schwarzschild[2][2] = 16.000000000000000;
    metric_schwarzschild[3][3] = 3.677581553054882;

    // christoffel symbols
    christoffel_symbols_schwarzschild[0][0][1] = 0.022058823529412;
    christoffel_symbols_schwarzschild[0][1][0] = 0.022058823529412;
    christoffel_symbols_schwarzschild[1][0][0] = 0.015937500000000;
    christoffel_symbols_schwarzschild[1][1][1] = -0.022058823529412;
    christoffel_symbols_schwarzschild[1][2][2] = -3.400000000000000;
    christoffel_symbols_schwarzschild[1][3][3] = -0.781486080024162;
    christoffel_symbols_schwarzschild[2][1][2] = 0.250000000000000;
    christoffel_symbols_schwarzschild[2][2][1] = 0.250000000000000;
    christoffel_symbols_schwarzschild[2][3][3] = -0.420735492403948;
    christoffel_symbols_schwarzschild[3][1][3] = 0.250000000000000;
    christoffel_symbols_schwarzschild[3][2][3] = 1.830487721712452;
    christoffel_symbols_schwarzschild[3][3][1] = 0.250000000000000;
    christoffel_symbols_schwarzschild[3][3][2] = 1.830487721712452;

    // riemanns tensor
    riemann_tensor_schwarzschild[0][1][0][1] = 0.011029411764706;
    riemann_tensor_schwarzschild[0][1][1][0] = -0.011029411764706;
    riemann_tensor_schwarzschild[0][2][0][2] = -0.075000000000000;
    riemann_tensor_schwarzschild[0][2][2][0] = 0.075000000000000;
    riemann_tensor_schwarzschild[0][3][0][3] = -0.017238663529945;
    riemann_tensor_schwarzschild[0][3][3][0] = 0.017238663529945;
    riemann_tensor_schwarzschild[1][0][0][1] = 0.007968750000000;
    riemann_tensor_schwarzschild[1][0][1][0] = -0.007968750000000;
    riemann_tensor_schwarzschild[1][2][1][2] = -0.075000000000000;
    riemann_tensor_schwarzschild[1][2][2][1] = 0.075000000000000;
    riemann_tensor_schwarzschild[1][3][1][3] = -0.017238663529945;
    riemann_tensor_schwarzschild[1][3][3][1] = 0.017238663529945;
    riemann_tensor_schwarzschild[2][0][0][2] = -0.003984375000000;
    riemann_tensor_schwarzschild[2][0][2][0] = 0.003984375000000;
    riemann_tensor_schwarzschild[2][1][1][2] = 0.005514705882353;
    riemann_tensor_schwarzschild[2][1][2][1] = -0.005514705882353;
    riemann_tensor_schwarzschild[2][3][2][3] = 0.034477327059890;
    riemann_tensor_schwarzschild[2][3][3][2] = -0.034477327059890;
    riemann_tensor_schwarzschild[3][0][0][3] = -0.003984375000000;
    riemann_tensor_schwarzschild[3][0][3][0] = 0.003984375000000;
    riemann_tensor_schwarzschild[3][1][1][3] = 0.005514705882353;
    riemann_tensor_schwarzschild[3][1][3][1] = -0.005514705882353;
    riemann_tensor_schwarzschild[3][2][2][3] = -0.150000000000000;
    riemann_tensor_schwarzschild[3][2][3][2] = 0.150000000000000;

    // ========== Prepare WeylSchwarzschild ========== 
    // coordinates
    y_weylschwarzschild[0] = 0.000000000000000;
    y_weylschwarzschild[1] = 0.000000000000000;
    y_weylschwarzschild[2] = 10.000000000000000;
    y_weylschwarzschild[3] = 2.000000000000000;

    // metric
    metric_weylschwarzschild[0][0] = -0.942869554762639;
    metric_weylschwarzschild[1][1] = 106.059209882086250;
    metric_weylschwarzschild[2][2] = 1.059710254638068;
    metric_weylschwarzschild[3][3] = 1.059710254638068;

    // christoffel symbols
    christoffel_symbols_weylschwarzschild[0][0][2] = 0.002827609949378;
    christoffel_symbols_weylschwarzschild[0][0][3] = 0.000565033002746;
    christoffel_symbols_weylschwarzschild[0][2][0] = 0.002827609949378;
    christoffel_symbols_weylschwarzschild[0][3][0] = 0.000565033002746;
    christoffel_symbols_weylschwarzschild[1][1][2] = 0.097172390050622;
    christoffel_symbols_weylschwarzschild[1][1][3] = -0.000565033002746;
    christoffel_symbols_weylschwarzschild[1][2][1] = 0.097172390050622;
    christoffel_symbols_weylschwarzschild[1][3][1] = -0.000565033002746;
    christoffel_symbols_weylschwarzschild[2][0][0] = 0.002515845555277;
    christoffel_symbols_weylschwarzschild[2][1][1] = -9.725325263218092;
    christoffel_symbols_weylschwarzschild[2][2][2] = -0.002750848792062;
    christoffel_symbols_weylschwarzschild[2][2][3] = -0.000533079143940;
    christoffel_symbols_weylschwarzschild[2][3][2] = -0.000533079143940;
    christoffel_symbols_weylschwarzschild[2][3][3] = 0.002750848792062;
    christoffel_symbols_weylschwarzschild[3][0][0] = 0.000502734038284;
    christoffel_symbols_weylschwarzschild[3][1][1] = 0.056550319831527;
    christoffel_symbols_weylschwarzschild[3][2][2] = 0.000533079143940;
    christoffel_symbols_weylschwarzschild[3][2][3] = -0.002750848792062;
    christoffel_symbols_weylschwarzschild[3][3][2] = -0.002750848792062;
    christoffel_symbols_weylschwarzschild[3][3][3] = -0.000533079143940;

    // ========== Prepare BachWeylRing ========== 
    // coordinates
    y_bachweylring[0] = 1.000000000000000;
    y_bachweylring[1] = 4.000000000000000;
    y_bachweylring[2] = 0.500000000000000;
    y_bachweylring[3] = 0.300000000000000;

    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}