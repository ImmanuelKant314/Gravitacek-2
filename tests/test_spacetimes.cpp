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

    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}