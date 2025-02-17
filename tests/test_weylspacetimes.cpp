#include <iostream>
#include <fstream>
#include <memory>

#include "gtest/gtest.h"
#include "gravitacek2/setup.hpp"
#include "gravitacek2/geomotion/spacetimes.hpp"

class WeylTestCase
{
public:
    std::shared_ptr<gr2::Weyl> spacetime;
    gr2::real y[4];
    gr2::real nu;
    gr2::real nu_rho, nu_z;
    gr2::real nu_rhorho, nu_rhoz, nu_zz;
    gr2::real lambda;
    gr2::real eps;

    std::string name;

    WeylTestCase(std::shared_ptr<gr2::Weyl> spacetime, std::string filename, std::string name, gr2::real eps):
    spacetime(spacetime), eps(eps), name(name)
    {
        std::ifstream file;
        file.open(filename);

        if (!file.is_open())
            throw std::runtime_error("Unable to open file " + filename);

        std::string line;
        int c;

        // position
        std::getline(file, line);
        c = sscanf(line.c_str(), "%Lf;%Lf;%Lf;%Lf%*s", y, y+1, y+2, y+3);
        if (c != 4)
            throw std::runtime_error("Unable to read position from file " + filename);

        // potential
        std::getline(file, line);
        c = sscanf(line.c_str(), "%Lf%*s", &nu);
        if (c != 1)
            throw std::runtime_error("Unable to read potential from file " + filename);

        // derivatives of potential 
        std::getline(file, line);
        c = sscanf(line.c_str(), "%Lf;%Lf%*s", &nu_rho, &nu_z);
        if (c != 2)
            throw std::runtime_error("Unable to read derivatives of potential from file " + filename);

        // second derivatives of potential
        std::getline(file, line);
        c = sscanf(line.c_str(), "%Lf;%Lf;%Lf%*s", &nu_rhorho, &nu_rhoz, &nu_zz);
        if (c != 3)
            throw std::runtime_error("Unable to read second derivatives of potential from file " + filename);

        // lambda
        std::getline(file, line);
        c = sscanf(line.c_str(), "%Lf%*s", &lambda);
        if (c != 1)
            throw std::runtime_error("Unable to read lambda from file " + filename);

        file.close();
    }
};

class GeneralWeylTest :
    public testing::TestWithParam<WeylTestCase> 
{

};

TEST_P(GeneralWeylTest, Potential)
{
    std::shared_ptr<gr2::Weyl> spacetime = GetParam().spacetime;
    gr2::real eps = GetParam().eps;
    gr2::real nu =  GetParam().nu;
    spacetime->calculate_nu(GetParam().y);
    EXPECT_NEAR(spacetime->get_nu(), nu, eps + eps*nu);
}

TEST_P(GeneralWeylTest, DerivativesOfPotential)
{
    std::shared_ptr<gr2::Weyl> spacetime = GetParam().spacetime;
    gr2::real eps = GetParam().eps;
    gr2::real nu_rho =  GetParam().nu_rho;
    gr2::real nu_z =  GetParam().nu_z;
    spacetime->calculate_nu1(GetParam().y);
    EXPECT_NEAR(spacetime->get_nu_rho(), nu_rho, eps + eps*nu_rho);
    EXPECT_NEAR(spacetime->get_nu_z(), nu_z, eps + eps*nu_z);
}

TEST_P(GeneralWeylTest, Lambda)
{
    std::shared_ptr<gr2::Weyl> spacetime = GetParam().spacetime;
    gr2::real eps = GetParam().eps;
    gr2::real lambda =  GetParam().lambda;
    spacetime->calculate_lambda_init(GetParam().y);
    EXPECT_NEAR(spacetime->get_lambda(), lambda, eps + eps*lambda);
}

void PrintTo(const WeylTestCase& testcase, std::ostream* os) {
    *os << "-";
}

std::string MyTestNameGenerator(const ::testing::TestParamInfo<WeylTestCase>& info) {
    return info.param.name;
}

auto test_cases = testing::Values(
    WeylTestCase(std::make_shared<gr2::WeylSchwarzschild>(0.3), "./test_weylspacetime/weylschwarzschild.txt", "WeylSchwarzschild", 1e-15),
    WeylTestCase(std::make_shared<gr2::BachWeylRing>(0.3, 5),  "./test_weylspacetime/bachweylring.txt", "BachWeylRing", 1e-15)
);

INSTANTIATE_TEST_SUITE_P(
    WeylTest,
    GeneralWeylTest,
    test_cases,
    MyTestNameGenerator
);

int main(int argc, char **argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}