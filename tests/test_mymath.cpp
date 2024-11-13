#include "gravitacek2/mymath.hpp"

#include "gtest/gtest.h"

TEST(elliptic_KE, K)
{
    int n = 10;
    gr2::real eps = 1e-15;

    gr2::real k[] = {0.000000000000000,
                    0.100000000000000,
                    0.200000000000000,
                    0.300000000000000,
                    0.400000000000000,
                    0.500000000000000,
                    0.600000000000000,
                    0.700000000000000,
                    0.800000000000000,
                    0.900000000000000};

    gr2::real K[] = {1.570796326794897,
                    1.574745561517356,
                    1.586867847454166,
                    1.608048619930513,
                    1.639999865864511,
                    1.685750354812596,
                    1.750753802915753,
                    1.845693998374724,
                    1.995302777664730,
                    2.280549138422770};
    
    gr2::real K_val, E_val;
    for (int i = 0; i < n; i++)
    {
        gr2::elliptic_KE(k[i], K_val, E_val, 1e-16);
        EXPECT_NEAR(K[i], K_val, eps);
    }
}

TEST(elliptic_KE, E)
{
    int n = 10;
    gr2::real eps = 1e-15;

    gr2::real k[] = {0.000000000000000,
                    0.100000000000000,
                    0.200000000000000,
                    0.300000000000000,
                    0.400000000000000,
                    0.500000000000000,
                    0.600000000000000,
                    0.700000000000000,
                    0.800000000000000,
                    0.900000000000000};

    gr2::real E[] = {1.570796326794897,
                    1.566861942021668,
                    1.554968546242529,
                    1.534833464923249,
                    1.505941612360040,
                    1.467462209339427,
                    1.418083394448724,
                    1.355661135571955,
                    1.276349943169906,
                    1.171697052781614};

    gr2::real K_val, E_val;
    for (int i = 0; i < n; i++)
    {
        gr2::elliptic_KE(k[i], K_val, E_val, 1e-16);
        EXPECT_NEAR(E[i], E_val, eps);
    }
}

TEST(romb5, IntegrateSinX)
{
    gr2::real eps=1e-13;
    EXPECT_NEAR(gr2::romb<5>(*sinl, 0, gr2::pi), 2, eps);
}

TEST(legendre_polynomials, Values)
{
    const int n = 6;
    gr2::real x = 0.27;
    gr2::real eps = 1e-14;

    gr2::real p_true[n] = {1.000000000000000,
                           0.270000000000000,
                          -0.390650000000000,
                          -0.355792500000000,
                           0.124875543750000,
                           0.345323514262500};
    gr2::real p_num[n] = {};

    gr2::legendre_polynomials(x, n, p_num);
    for (int i = 0; i < n; i++)
    {
        EXPECT_NEAR(p_true[i], p_num[i], eps);
    }
}

TEST(legendre_polynomials1, Values)
{
    const int n = 6;
    gr2::real x = 0.27;
    gr2::real eps = 1e-14;

    gr2::real p_true[n] = {1.000000000000000,
                           0.270000000000000,
                          -0.390650000000000,
                          -0.355792500000000,
                           0.124875543750000,
                           0.345323514262500};
    gr2::real p_num[n] = {};
    gr2::real p_extra[n];

    gr2::legendre_polynomials1(x, n, p_num, p_extra);
    for (int i = 0; i < n; i++)
    {
        EXPECT_NEAR(p_true[i], p_num[i], eps);
    }
}

TEST(legendre_polynomials1, Derivatives)
{
    const int n = 6;
    gr2::real x = 0.27;
    gr2::real eps = 1e-14;

    gr2::real p_true[n] = {0.000000000000000,
                           1.000000000000000,
                           0.810000000000000,
                          -0.953250000000000,
                          -1.680547500000000,
                           0.170629893750000};
    gr2::real p_num[n] = {};
    gr2::real p_extra[n];

    gr2::legendre_polynomials1(x, n, p_extra, p_num);
    for (int i = 0; i < n; i++)
    {
        EXPECT_NEAR(p_true[i], p_num[i], eps);
    }
}

TEST(special_function_Q2n, Values)
{
    const int n = 6;
    gr2::real x = 0.27;
    gr2::real eps = 1e-14;

    gr2::real q_true[n] = {1.307084492332630,
                          -0.391471935402888,
                           0.173008207890096,
                          -0.084723372649464,
                           0.043524450781026,
                          -0.022988984539059};
    gr2::real q_num[n];

    gr2::special_function_Q2n(x, n, q_num);
    for (int i = 0; i < n; i++)
    {
        EXPECT_NEAR(q_true[i], q_num[i], eps);
    }
}

TEST(special_function_Q2n1, Values)
{
    const int n = 6;
    gr2::real x = 0.27;
    gr2::real eps = 1e-14;

    gr2::real q_true[n] = {1.307084492332630,
                          -0.391471935402888,
                           0.173008207890096,
                          -0.084723372649464,
                           0.043524450781026,
                          -0.022988984539059};
    gr2::real q_num[n];
    gr2::real q_help[n];

    gr2::special_function_Q2n1(x, n, q_num, q_help);
    for (int i = 0; i < n; i++)
    {
        EXPECT_NEAR(q_true[i], q_num[i], eps);
    }
}

TEST(special_function_Q2n, Derivatives)
{
    const int n = 6;
    gr2::real x = 0.27;
    gr2::real eps = 1e-14;

    gr2::real q_true[n] = {-0.932053313449529,
                            1.009208247761040,
                           -0.777395362047414,
                            0.543718725188478,
                           -0.363199801841301,
                            0.236171352375938};
    gr2::real q_num[n];
    gr2::real q_help[n];

    gr2::special_function_Q2n1(x, n, q_help, q_num);
    for (int i = 0; i < n; i++)
    {
        EXPECT_NEAR(q_true[i], q_num[i], eps);
    }
}

int main(int argc, char **argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}