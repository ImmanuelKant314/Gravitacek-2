#include <iostream>
#include <fstream>

#include "gtest/gtest.h"
#include "gravitacek2/setup.hpp"
#include "gravitacek2/geomotion/spacetimes.hpp"

// ========== prepare Schwarzschild ========== 

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

gr2::real y_weylschwarzschild[4] = {};
gr2::real metric_weylschwarzschild[4][4] = {};
gr2::real christoffel_symbols_weylschwarzschild[4][4][4] = {};
gr2::real riemann_tensor_weylschwarzschild[4][4][4][4] = {};

TEST(WeylSchwarzschild, potencial)
{
    gr2::WeylSchwarzschild weylschwarzschild = gr2::WeylSchwarzschild(0.3);
    weylschwarzschild.calculate_nu(y_weylschwarzschild);
    EXPECT_NEAR(weylschwarzschild.get_nu(), -0.029413667983693322, eps);
}

TEST(WeylSchwarzchild, DerivativesOfPotential)
{
    gr2::WeylSchwarzschild weylschwarzschild = gr2::WeylSchwarzschild(0.3);
    weylschwarzschild.calculate_nu1(y_weylschwarzschild);
    EXPECT_NEAR(weylschwarzschild.get_nu(), -0.029413667983693322, eps);
    EXPECT_NEAR(weylschwarzschild.get_nu_rho(), 0.0028276099493778893, eps);
    EXPECT_NEAR(weylschwarzschild.get_nu_z(), 0.0005650330027458911, eps);
}

TEST(WeylSchwarzschild, Lambda)
{
    gr2::WeylSchwarzschild weylschwarzschild = gr2::WeylSchwarzschild(0.3);
    weylschwarzschild.calculate_lambda(y_weylschwarzschild);
    EXPECT_NEAR(weylschwarzschild.get_lambda(), -0.000415904944973, eps);
}

TEST(WeylSchwarzschild, Metric)
{
    gr2::WeylSchwarzschild weylschwarzschild = gr2::WeylSchwarzschild(0.3);
    weylschwarzschild.calculate_metric(y_weylschwarzschild);
    for (int i = 0; i < 4; i++)
        for (int j = 0; j < 4; j++)
            EXPECT_NEAR(weylschwarzschild.get_metric()[i][j], metric_weylschwarzschild[i][j], eps);
}

TEST(WeylSchwarzschild, ChristoffelSymbols)
{
    gr2::WeylSchwarzschild weylschwarzschild = gr2::WeylSchwarzschild(0.3);
    weylschwarzschild.calculate_christoffel_symbols(y_weylschwarzschild);
    for (int i = 0; i < 4; i++)
        for (int j = 0; j < 4; j++)
            for(int k = 0; k < 4; k++)
                EXPECT_NEAR(weylschwarzschild.get_christoffel_symbols()[i][j][k], christoffel_symbols_weylschwarzschild[i][j][k], 4*eps);
}

//TODO: WeylSchwarzschild - Riemannův tensor

int main(int argc, char **argv)
{
    // ========== prepare schwarzschild ========== 
    std::ifstream file("data/schwarzschild.txt");

    // vector
    while (getline(file,line))
    {
        if (line == "")
            break;
        if (line[0] == '#')
            continue;

        int i = std::stoi(line.substr(0, 1));
        gr2::real value = std::stold(line.substr(2));
        y_schwarzschild[i] = value;
    }

    // metric
    while (getline(file,line))
    {
        if (line == "")
            break;
        if (line[0] == '#')
            continue;

        int i = std::stoi(line.substr(0, 1));
        int j = std::stoi(line.substr(2, 1));
        gr2::real value = std::stold(line.substr(4));
        metric_schwarzschild[i][j] = value;
    }

    // christoffel symbols
    while (getline(file,line))
    {
        if (line == "")
            break;
        if (line[0] == '#')
            continue;

        int i = std::stoi(line.substr(0, 1));
        int j = std::stoi(line.substr(2, 1));
        int k = std::stoi(line.substr(4, 1));
        gr2::real value = std::stold(line.substr(6));
        christoffel_symbols_schwarzschild[i][j][k] = value;
    }

    // riemanns tensor
    while (getline(file,line))
    {
        if (line == "")
            break;
        if (line[0] == '#')
            continue;

        int i = std::stoi(line.substr(0, 1));
        int j = std::stoi(line.substr(2, 1));
        int k = std::stoi(line.substr(4, 1));
        int l = std::stoi(line.substr(6, 1));
        gr2::real value = std::stold(line.substr(8));
        riemann_tensor_schwarzschild[i][j][k][l] = value;
    }
    
    file.close();

    y_weylschwarzschild[0] = 1.000000000000000;
    y_weylschwarzschild[1] = 3.000000000000000;
    y_weylschwarzschild[2] = 10.000000000000000;
    y_weylschwarzschild[3] = 2.000000000000000;

    metric_weylschwarzschild[0][0] = -0.942869554762639;
    metric_weylschwarzschild[1][1] = 106.059209882086250;
    metric_weylschwarzschild[2][2] = 1.059710254638068;
    metric_weylschwarzschild[3][3] = 1.059710254638068;

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

    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}