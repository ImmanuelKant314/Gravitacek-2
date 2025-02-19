#include <iostream>
#include <fstream>
#include <memory>

#include "gtest/gtest.h"
#include "gravitacek2/setup.hpp"
#include "gravitacek2/geomotion/spacetimes.hpp"

class MPTestCase
{
public:
    std::shared_ptr<gr2::MajumdarPapapetrouWeyl> spacetime;
    gr2::real y[4];
    gr2::real N;
    gr2::real N_inv;
    gr2::real N_rho, N_z;
    gr2::real N_rhorho, N_rhoz, N_zz;
    gr2::real eps;

    std::string name;

    MPTestCase(std::shared_ptr<gr2::MajumdarPapapetrouWeyl> spacetime, std::string filename, std::string name, gr2::real eps):
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

        // lapse
        std::getline(file, line);
        c = sscanf(line.c_str(), "%Lf%*s", &N);
        if (c != 1)
            throw std::runtime_error("Unable to read lapse function from file " + filename);

        // inverse of lapse
        std::getline(file, line);
        c = sscanf(line.c_str(), "%Lf%*s", &N_inv);
        if (c != 1)
            throw std::runtime_error("Unable to read inverse of lapse function from file " + filename);

        // derivatives of lapse
        std::getline(file, line);
        c = sscanf(line.c_str(), "%Lf;%Lf%*s", &N_rho, &N_z);
        if (c != 2)
            throw std::runtime_error("Unable to read derivatives of lapse from file " + filename);

        // second derivatives of lapse
        std::getline(file, line);
        c = sscanf(line.c_str(), "%Lf;%Lf;%Lf%*s", &N_rhorho, &N_rhoz, &N_zz);
        if (c != 3)
            throw std::runtime_error("Unable to read second derivatives of lapse from file " + filename);

        file.close();
    }
};

class GeneralMPTest :
    public testing::TestWithParam<MPTestCase> 
{

};

TEST_P(GeneralMPTest, Lapse)
{
    std::shared_ptr<gr2::MajumdarPapapetrouWeyl> spacetime = GetParam().spacetime;
    gr2::real eps = GetParam().eps;
    gr2::real N =  GetParam().N;
    spacetime->calculate_N(GetParam().y);
    EXPECT_NEAR(spacetime->get_N(), N, eps + eps*std::llabs(N));
}

TEST_P(GeneralMPTest, InvLapse)
{
    std::shared_ptr<gr2::MajumdarPapapetrouWeyl> spacetime = GetParam().spacetime;
    gr2::real eps = GetParam().eps;
    gr2::real N_inv =  GetParam().N_inv;
    spacetime->calculate_N(GetParam().y);
    EXPECT_NEAR(spacetime->get_invN(), N_inv, eps + eps*std::llabs(N_inv));
}

TEST_P(GeneralMPTest, DerivativesOfLapse)
{
    std::shared_ptr<gr2::MajumdarPapapetrouWeyl> spacetime = GetParam().spacetime;
    gr2::real eps = GetParam().eps;
    gr2::real N =  GetParam().N;
    gr2::real N_rho =  GetParam().N_rho;
    gr2::real N_z =  GetParam().N_z;
    spacetime->calculate_N1(GetParam().y);
    EXPECT_NEAR(spacetime->get_N(), N, eps + eps*std::llabs(N));
    EXPECT_NEAR(spacetime->get_N_rho(), N_rho, eps + eps*std::llabs(N_rho));
    EXPECT_NEAR(spacetime->get_N_z(), N_z, eps + eps*std::llabs(N_z));
}

TEST_P(GeneralMPTest, SecondDerivativesOfLapse)
{
    std::shared_ptr<gr2::MajumdarPapapetrouWeyl> spacetime = GetParam().spacetime;
    gr2::real eps = GetParam().eps;
    gr2::real N = GetParam().N;
    gr2::real N_rho = GetParam().N_rho;
    gr2::real N_z = GetParam().N_z;
    gr2::real N_rhorho = GetParam().N_rhorho;
    gr2::real N_rhoz = GetParam().N_rhoz;
    gr2::real N_zz = GetParam().N_zz;
    spacetime->calculate_N2(GetParam().y);
    EXPECT_NEAR(spacetime->get_N(), N, eps + eps*std::llabs(N));
    EXPECT_NEAR(spacetime->get_N_rho(), N_rho, eps + eps*std::llabs(N_rho));
    EXPECT_NEAR(spacetime->get_N_z(), N_z, eps + eps*std::llabs(N_z));
    EXPECT_NEAR(spacetime->get_N_rhorho(), N_rhorho, eps + eps*std::llabs(N_rhorho));
    EXPECT_NEAR(spacetime->get_N_rhoz(), N_rhoz, eps + eps*std::llabs(N_rhoz));
    EXPECT_NEAR(spacetime->get_N_zz(), N_zz, eps + eps*std::llabs(N_zz));
}

void PrintTo(const MPTestCase& testcase, std::ostream* os) {
    *os << "0";
}

std::string MyTestNameGenerator(const ::testing::TestParamInfo<MPTestCase>& info) {
    return info.param.name;
}

std::string folder = "./test_majumdarpapapetrouwspacetime/";

auto test_cases = testing::Values(
    MPTestCase(std::make_shared<gr2::ReissnerNordstromMPW>(0.3), folder + "reissnernordstrom.txt", "ReissnerNordstrom", 1e-12)
);

INSTANTIATE_TEST_SUITE_P(
    MPTest,
    GeneralMPTest,
    test_cases,
    MyTestNameGenerator
);

int main(int argc, char **argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}