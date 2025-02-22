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
    gr2::real N_inv;
    gr2::real N_inv_rho, N_inv_z;
    gr2::real N_inv_rhorho, N_inv_rhoz, N_inv_zz;
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

        // inverse of lapse
        std::getline(file, line);
        c = sscanf(line.c_str(), "%Lf%*s", &N_inv);
        if (c != 1)
            throw std::runtime_error("Unable to read inverse of lapse function from file " + filename);

        // derivatives of inverse lapse
        std::getline(file, line);
        c = sscanf(line.c_str(), "%Lf;%Lf%*s", &N_inv_rho, &N_inv_z);
        if (c != 2)
            throw std::runtime_error("Unable to read derivatives of inverse lapse from file " + filename);

        // second derivatives of inverse lapse
        std::getline(file, line);
        c = sscanf(line.c_str(), "%Lf;%Lf;%Lf%*s", &N_inv_rhorho, &N_inv_rhoz, &N_inv_zz);
        if (c != 3)
            throw std::runtime_error("Unable to read second derivatives of lapse inverse from file " + filename);

        file.close();
    }
};

class GeneralMPTest :
    public testing::TestWithParam<MPTestCase> 
{

};

TEST_P(GeneralMPTest, InverseLapse)
{
    std::shared_ptr<gr2::MajumdarPapapetrouWeyl> spacetime = GetParam().spacetime;
    gr2::real eps = GetParam().eps;
    gr2::real N_inv =  GetParam().N_inv;
    spacetime->calculate_N_inv(GetParam().y);
    EXPECT_NEAR(spacetime->get_N_inv(), N_inv, eps + eps*std::llabs(N_inv));
}

TEST_P(GeneralMPTest, DerivativesOfInverseLapse)
{
    std::shared_ptr<gr2::MajumdarPapapetrouWeyl> spacetime = GetParam().spacetime;
    gr2::real eps = GetParam().eps;
    gr2::real N_inv =  GetParam().N_inv;
    gr2::real N_inv_rho =  GetParam().N_inv_rho;
    gr2::real N_inv_z =  GetParam().N_inv_z;
    spacetime->calculate_N_inv1(GetParam().y);
    EXPECT_NEAR(spacetime->get_N_inv(), N_inv, eps + eps*std::llabs(N_inv));
    EXPECT_NEAR(spacetime->get_N_inv_rho(), N_inv_rho, eps + eps*std::llabs(N_inv_rho));
    EXPECT_NEAR(spacetime->get_N_inv_z(), N_inv_z, eps + eps*std::llabs(N_inv_z));
}

TEST_P(GeneralMPTest, SecondDerivativesOfLapse)
{
    std::shared_ptr<gr2::MajumdarPapapetrouWeyl> spacetime = GetParam().spacetime;
    gr2::real eps = GetParam().eps;
    gr2::real N_inv = GetParam().N_inv;
    gr2::real N_inv_rho = GetParam().N_inv_rho;
    gr2::real N_inv_z = GetParam().N_inv_z;
    gr2::real N_inv_rhorho = GetParam().N_inv_rhorho;
    gr2::real N_inv_rhoz = GetParam().N_inv_rhoz;
    gr2::real N_inv_zz = GetParam().N_inv_zz;
    spacetime->calculate_N_inv2(GetParam().y);
    EXPECT_NEAR(spacetime->get_N_inv(), N_inv, eps + eps*std::llabs(N_inv));
    EXPECT_NEAR(spacetime->get_N_inv_rho(), N_inv_rho, eps + eps*std::llabs(N_inv_rho));
    EXPECT_NEAR(spacetime->get_N_inv_z(), N_inv_z, eps + eps*std::llabs(N_inv_z));
    EXPECT_NEAR(spacetime->get_N_inv_rhorho(), N_inv_rhorho, eps + eps*std::llabs(N_inv_rhorho));
    EXPECT_NEAR(spacetime->get_N_inv_rhoz(), N_inv_rhoz, eps + eps*std::llabs(N_inv_rhoz));
    EXPECT_NEAR(spacetime->get_N_inv_zz(), N_inv_zz, eps + eps*std::llabs(N_inv_zz));
}

void PrintTo(const MPTestCase& testcase, std::ostream* os) {
    *os << "0";
}

std::string MyTestNameGenerator(const ::testing::TestParamInfo<MPTestCase>& info) {
    return info.param.name;
}

std::string folder = "./test_majumdarpapapetrouwspacetime/";

std::vector<std::shared_ptr<gr2::MajumdarPapapetrouWeyl>> vec_rn_mp = {std::make_shared<gr2::ReissnerNordstromMPW>(1.0), std::make_shared<gr2::MajumdarPapapetrouRing>(1.0, 5.0)};
std::shared_ptr<gr2::CombinedMPW> rnmp = std::make_shared<gr2::CombinedMPW>(vec_rn_mp);

auto test_cases = testing::Values(
    MPTestCase(std::make_shared<gr2::ReissnerNordstromMPW>(0.3), folder + "reissnernordstrom.txt", "ReissnerNordstrom", 1e-12),
    MPTestCase(std::make_shared<gr2::MajumdarPapapetrouRing>(0.3, 5), folder + "majumdarpapapetrouring.txt", "MajumdarPapapetrouRing", 1e-12),
    MPTestCase(rnmp, folder + "rnmpr.txt", "ReissnerNordstromMajumdarPapapetrouRing", 1e-12)
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