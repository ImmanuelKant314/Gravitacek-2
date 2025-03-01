#include <iostream>
#include <fstream>

#include "gtest/gtest.h"
#include "gravitacek2/setup.hpp"
#include "gravitacek2/geomotion/geomotion.hpp"
#include "gravitacek2/geomotion/spacetimes.hpp"

class SpacetimeTestCase
{
public:
    std::shared_ptr<gr2::GeoMotion> spacetime;
    gr2::real y[4];
    gr2::real metric[4][4];
    gr2::real christoffel_symbols[4][4][4];
    gr2::real riemann_tensor[4][4][4][4];

    gr2::real eps;

    std::string name;
    SpacetimeTestCase(std::shared_ptr<gr2::GeoMotion> spacetime, std::string filename, std::string name, gr2::real eps):
    spacetime(spacetime), eps(eps), name(name), metric(), christoffel_symbols(), riemann_tensor()
    {
        std::ifstream file;
        file.open(filename);

        if (!file.is_open())
            throw std::runtime_error("Unable to open file " + filename);

        std::string line;
        int c;
        int i, j, k, l;
        gr2::real value;

        // position
        std::getline(file, line);
        c = sscanf(line.c_str(), "%Lf;%Lf;%Lf;%Lf%*s", y, y+1, y+2, y+3);
        if (c != 4)
            throw std::runtime_error("Unable to read position from file " + filename);

        std::getline(file, line); // read line containing =====

        // metric
        while(true)
        {
            std::getline(file, line);
            if (line[0] == '=')
                break;
            if (file.eof())
                throw std::runtime_error("File " + filename + "does not contain all information for the test");
            c = sscanf(line.c_str(), "%d;%d;%Lf%*s", &i, &j, &value);
            if (c != 3)
                throw std::runtime_error("Unable to read metric component from file " + filename);
            this->metric[i][j] = value;
        }

        // Christoffel symbols
        while(true)
        {
            std::getline(file, line);
            if (line[0] == '=')
                break;
            if (file.eof())
                throw std::runtime_error("File " + filename + "does not contain all information for the test");
            c = sscanf(line.c_str(), "%d;%d;%d;%Lf%*s", &i, &j, &k, &value);
            if (c != 4)
                throw std::runtime_error("Unable to read Christoffel symbol component from file " + filename);
            this->christoffel_symbols[i][j][k] = value;
        }

        // Riemann tensor
        while(true)
        {
            std::getline(file, line);
            if (line[0] == '=')
                break;
            c = sscanf(line.c_str(), "%d;%d;%d;%d;%Lf%*s", &i, &j, &k, &l, &value);
            if (c != 5)
                throw std::runtime_error("Unable to read Christoffel symbol component from file " + filename);
            this->riemann_tensor[i][j][k][l] = value;
            if (file.eof())
                break;

        }
        file.close();
    }
};

class GeneralSpacetimeTest :
    public testing::TestWithParam<SpacetimeTestCase>
{

};

TEST_P(GeneralSpacetimeTest, Metric)
{
    std::shared_ptr<gr2::GeoMotion> spacetime = GetParam().spacetime;
    gr2::real eps = GetParam().eps;
    spacetime->calculate_metric(GetParam().y);
    for (int i = 0; i <4; i++)
        for (int j = 0; j < 4; j++)
        {
            gr2::real value = GetParam().metric[i][j];
            EXPECT_NEAR(spacetime->get_metric()[i][j], value, eps + std::abs(value)*eps);
        }
}

TEST_P(GeneralSpacetimeTest, ChristoffelSymbols)
{
    std::shared_ptr<gr2::GeoMotion> spacetime = GetParam().spacetime;
    gr2::real eps = GetParam().eps;
    spacetime->calculate_christoffel_symbols(GetParam().y);
    for (int i = 0; i <4; i++)
        for (int j = 0; j < 4; j++)
            for (int k = 0; k < 4; k++)
            {
                gr2::real value = GetParam().christoffel_symbols[i][j][k];
                EXPECT_NEAR(spacetime->get_christoffel_symbols()[i][j][k], value, eps + std::abs(value)*eps);
            }
}

TEST_P(GeneralSpacetimeTest, RiemannTensor)
{
    std::shared_ptr<gr2::GeoMotion> spacetime = GetParam().spacetime;
    gr2::real eps = GetParam().eps;
    spacetime->calculate_riemann_tensor(GetParam().y);
    for (int i = 0; i <4; i++)
        for (int j = 0; j < 4; j++)
            for (int k = 0; k < 4; k++)
                for (int l = 0; l < 4; l++)
                {
                    gr2::real value = GetParam().riemann_tensor[i][j][k][l];
                    EXPECT_NEAR(spacetime->get_riemann_tensor()[i][j][k][l], value, eps + std::abs(value)*eps);
                }
}

void PrintTo(const SpacetimeTestCase& testcase, std::ostream* os) {
    *os << "0";
}

std::string MyTestNameGenerator(const ::testing::TestParamInfo<SpacetimeTestCase>& info) {
    return info.param.name;
}

std::string folder = "./test_spacetime/";

auto test_cases = testing::Values(
    SpacetimeTestCase(std::make_shared<gr2::Schwarzschild>(0.3), folder + "schwarzschild.txt", "Schwarzschild", 1e-14),
    SpacetimeTestCase(std::make_shared<gr2::WeylSchwarzschild>(0.3, gr2::exact, gr2::exact), folder + "weyl.txt", "Weyl", 1e-14),
    SpacetimeTestCase(std::make_shared<gr2::ReissnerNordstromMPW>(0.3), folder + "mp.txt", "MajumdarPapapetrouWeyl", 1e-14)
);

INSTANTIATE_TEST_SUITE_P(
    SpacetimeTest,
    GeneralSpacetimeTest,
    test_cases,
    MyTestNameGenerator
);

int main(int argc, char **argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}