#include "gravitacek2/integrator/steppers.hpp"

#include <cmath>

namespace gr2
{
    DoPr853::DoPr853() : StepperBase(), k1(nullptr), k2(nullptr), k3(nullptr), k4(nullptr), k5(nullptr), k6(nullptr), k7(nullptr), k8(nullptr), k9(nullptr), k10(nullptr), k11(nullptr), k12(nullptr)
    {}

    DoPr853::~DoPr853()
    {
        delete[] k1, k2, k3, k4, k5, k6, k7, k8, k9, k10, k11, k12;
        delete[] k_help;
    }

    void DoPr853::set_OdeSystem(OdeSystem& ode)
    {
        int old_n = n;
        this->StepperBase::set_OdeSystem(ode);
        if(old_n != n)
        {
            delete[] k1, k2, k3, k4, k5, k6, k7, k8, k9, k10, k11, k12;
            k1 = new real[n];
            k2 = new real[n];
            k3 = new real[n];
            k4 = new real[n];
            k5 = new real[n];
            k6 = new real[n];
            k7 = new real[n];
            k8 = new real[n];
            k9 = new real[n];
            k10 = new real[n];
            k11 = new real[n];
            k12 = new real[n];
            k_help = new real[n];
        }
    }

    void DoPr853::reset()
    {
    }

    void DoPr853::step(const real &t, real y[], const real &h, const real dydt_in[], real dydt_out[]) 
    {
        int i = 0;
        static const real c2 = 0.526001519587677318785587544488e-01;
        static const real c3 = 0.789002279381515978178381316732e-01;
        static const real c4 = 0.118350341907227396726757197510e+00;
        static const real c5 = 0.281649658092772603273242802490e+00;
        static const real c6 = 0.333333333333333333333333333333e+00;
        static const real c7 = 0.25e+00;
        static const real c8 = 0.307692307692307692307692307692e+00;
        static const real c9 = 0.651282051282051282051282051282e+00;
        static const real c10 = 0.6e+00;
        static const real c11 = 0.857142857142857142857142857142e+00;
        static const real c14 = 0.1e+00;
        static const real c15 = 0.2e+00;
        static const real c16 = 0.777777777777777777777777777778e+00;

        static const real b1 = 5.42937341165687622380535766363e-2;
        static const real b6 = 4.45031289275240888144113950566e0;
        static const real b7 = 1.89151789931450038304281599044e0;
        static const real b8 = -5.8012039600105847814672114227e0;
        static const real b9 = 3.1116436695781989440891606237e-1;
        static const real b10 = -1.52160949662516078556178806805e-1;
        static const real b11 = 2.01365400804030348374776537501e-1;
        static const real b12 = 4.47106157277725905176885569043e-2;
        static const real bhh1 = 0.244094488188976377952755905512e+00;
        static const real bhh2 = 0.733846688281611857341361741547e+00;
        static const real bhh3 = 0.220588235294117647058823529412e-01;

        static const real er1 = 0.1312004499419488073250102996e-01;
        static const real er6 = -0.1225156446376204440720569753e+01;
        static const real er7 = -0.4957589496572501915214079952e+00;
        static const real er8 = 0.1664377182454986536961530415e+01;
        static const real er9 = -0.3503288487499736816886487290e+00;
        static const real er10 = 0.3341791187130174790297318841e+00;
        static const real er11 = 0.8192320648511571246570742613e-01;
        static const real er12 = -0.2235530786388629525884427845e-01;

        static const real a21 = 5.26001519587677318785587544488e-2;

        static const real a31 = 1.97250569845378994544595329183e-2;
        static const real a32 = 5.91751709536136983633785987549e-2;

        static const real a41 = 2.95875854768068491816892993775e-2;
        static const real a43 = 8.87627564304205475450678981324e-2;

        static const real a51 = 2.41365134159266685502369798665e-1;
        static const real a53 = -8.84549479328286085344864962717e-1;
        static const real a54 = 9.24834003261792003115737966543e-1;

        static const real a61 = 3.7037037037037037037037037037e-2;
        static const real a64 = 1.70828608729473871279604482173e-1;
        static const real a65 = 1.25467687566822425016691814123e-1;

        static const real a71 = 3.7109375e-2;
        static const real a74 = 1.70252211019544039314978060272e-1;
        static const real a75 = 6.02165389804559606850219397283e-2;
        static const real a76 = -1.7578125e-2;

        static const real a81 = 3.70920001185047927108779319836e-2;
        static const real a84 = 1.70383925712239993810214054705e-1;
        static const real a85 = 1.07262030446373284651809199168e-1;
        static const real a86 = -1.53194377486244017527936158236e-2;
        static const real a87 = 8.27378916381402288758473766002e-3;

        static const real a91 = 6.24110958716075717114429577812e-1;
        static const real a94 = -3.36089262944694129406857109825e0;
        static const real a95 = -8.68219346841726006818189891453e-1;
        static const real a96 = 2.75920996994467083049415600797e1;
        static const real a97 = 2.01540675504778934086186788979e1;
        static const real a98 = -4.34898841810699588477366255144e1;

        static const real a101 = 4.77662536438264365890433908527e-1;
        static const real a104 = -2.48811461997166764192642586468e0;
        static const real a105 = -5.90290826836842996371446475743e-1;
        static const real a106 = 2.12300514481811942347288949897e1;
        static const real a107 = 1.52792336328824235832596922938e1;
        static const real a108 = -3.32882109689848629194453265587e1;
        static const real a109 = -2.03312017085086261358222928593e-2;

        static const real a111 = -9.3714243008598732571704021658e-1;
        static const real a114 = 5.18637242884406370830023853209e0;
        static const real a115 = 1.09143734899672957818500254654e0;
        static const real a116 = -8.14978701074692612513997267357e0;
        static const real a117 = -1.85200656599969598641566180701e1;
        static const real a118 = 2.27394870993505042818970056734e1;
        static const real a119 = 2.49360555267965238987089396762e0;
        static const real a1110 = -3.0467644718982195003823669022e0;

        static const real a121 = 2.27331014751653820792359768449e0;
        static const real a124 = -1.05344954667372501984066689879e1;
        static const real a125 = -2.00087205822486249909675718444e0;
        static const real a126 = -1.79589318631187989172765950534e1;
        static const real a127 = 2.79488845294199600508499808837e1;
        static const real a128 = -2.85899827713502369474065508674e0;
        static const real a129 = -8.87285693353062954433549289258e0;
        static const real a1210 = 1.23605671757943030647266201528e1;
        static const real a1211 = 6.43392746015763530355970484046e-1;

        if (dydt_in)
        {
            // copy dydt
            for (i = 0; i < n; i++)
                k1[i] = dydt_in[i];
        }
        else
        {
            // first correction
            for (i = 0; i < n; i++)
                yt[i] = y[i];
            ode->function(t, yt, k1);
        }

        // 2. corection
        for (i = 0; i < n; i++)
            yt[i] = y[i] + h * (a21 * k1[i]);
        ode->function(t + c2 * h, yt, k2);

        // 3. corection
        for (i = 0; i < n; i++)
            yt[i] = y[i] + h * (a31 * k1[i] + a32 * k2[i]);
        ode->function(t + c3 * h, yt, k3);

        // 4. corection
        for (i = 0; i < n; i++)
            yt[i] = y[i] + h * (a41 * k1[i] + a43 * k3[i]);
        ode->function(t + c4 * h, yt, k4);

        // 5. corection
        for (i = 0; i < n; i++)
            yt[i] = y[i] + h * (a51 * k1[i] + a53 * k3[i] + a54 * k4[i]);
        ode->function(t + c5 * h, yt, k5);

        // 6. corection
        for (i = 0; i < n; i++)
            yt[i] = y[i] + h * (a61 * k1[i] + a64 * k4[i] + a65 * k5[i]);
        ode->function(t + c6 * h, yt, k6);

        // 7. corection
        for (i = 0; i < n; i++)
            yt[i] = y[i] + h * (a71 * k1[i] + a74 * k4[i] + a75 * k5[i] + a76 * k6[i]);
        ode->function(t + c7 * h, yt, k7);

        // 8. corection
        for (i = 0; i < n; i++)
            yt[i] = y[i] + h * (a81 * k1[i] + a84 * k4[i] + a85 * k5[i] + a86 * k6[i] + a87 * k7[i]);
        ode->function(t + c8 * h, yt, k8);

        // 9. corection
        for (i = 0; i < n; i++)
            yt[i] = y[i] + h * (a91 * k1[i] + a94 * k4[i] + a95 * k5[i] + a96 * k6[i] + a97 * k7[i] + a98 * k8[i]);
        ode->function(t + c9 * h, yt, k9);

        // 10. corection
        for (i = 0; i < n; i++)
            yt[i] = y[i] + h * (a101 * k1[i] + a104 * k4[i] + a105 * k5[i] + a106 * k6[i] + a107 * k7[i] + a108 * k8[i] + a109 * k9[i]);
        ode->function(t + c10 * h, yt, k10);

        // 11. corection
        for (i = 0; i < n; i++)
            yt[i] = y[i] + h * (a111 * k1[i] + a114 * k4[i] + a115 * k5[i] + a116 * k6[i] + a117 * k7[i] + a118 * k8[i] + a119 * k9[i] + a1110 * k10[i]);
        ode->function(t + c11 * h, yt, k11);

        // 12. corection
        for (i = 0; i < n; i++)
            yt[i] = y[i] + h * (a121 * k1[i] + a124 * k4[i] + a125 * k5[i] + a126 * k6[i] + a127 * k7[i] + a128 * k8[i] + a129 * k9[i] + a1210 * k10[i] + a1211 * k11[i]);
        ode->function(t, yt, k12);

        // final value
        for (i = 0; i < n; i++)
        {
            k_help[i] = (b1 * k1[i] + b6 * k6[i] + b7 * k7[i] + b8 * k8[i] + b9 * k9[i] + b10 * k10[i] + b11 * k11[i] + b12 * k12[i]);
            y[i] = y[i] + h * k_help[i];
        }

        if (dydt_out)
        {
            ode->function(t+h, y, dydt_out);
        }
    }

    void DoPr853::step_err(const real &t, real y[], const real &h, real err[], const real dydt_in[], real dydt_out[])
    {
        int i = 0;
        static const real c2 = 0.526001519587677318785587544488e-01;
        static const real c3 = 0.789002279381515978178381316732e-01;
        static const real c4 = 0.118350341907227396726757197510e+00;
        static const real c5 = 0.281649658092772603273242802490e+00;
        static const real c6 = 0.333333333333333333333333333333e+00;
        static const real c7 = 0.25e+00;
        static const real c8 = 0.307692307692307692307692307692e+00;
        static const real c9 = 0.651282051282051282051282051282e+00;
        static const real c10 = 0.6e+00;
        static const real c11 = 0.857142857142857142857142857142e+00;
        static const real c14 = 0.1e+00;
        static const real c15 = 0.2e+00;
        static const real c16 = 0.777777777777777777777777777778e+00;

        static const real b1 = 5.42937341165687622380535766363e-2;
        static const real b6 = 4.45031289275240888144113950566e0;
        static const real b7 = 1.89151789931450038304281599044e0;
        static const real b8 = -5.8012039600105847814672114227e0;
        static const real b9 = 3.1116436695781989440891606237e-1;
        static const real b10 = -1.52160949662516078556178806805e-1;
        static const real b11 = 2.01365400804030348374776537501e-1;
        static const real b12 = 4.47106157277725905176885569043e-2;
        static const real bhh1 = 0.244094488188976377952755905512e+00;
        static const real bhh2 = 0.733846688281611857341361741547e+00;
        static const real bhh3 = 0.220588235294117647058823529412e-01;

        static const real er1 = 0.1312004499419488073250102996e-01;
        static const real er6 = -0.1225156446376204440720569753e+01;
        static const real er7 = -0.4957589496572501915214079952e+00;
        static const real er8 = 0.1664377182454986536961530415e+01;
        static const real er9 = -0.3503288487499736816886487290e+00;
        static const real er10 = 0.3341791187130174790297318841e+00;
        static const real er11 = 0.8192320648511571246570742613e-01;
        static const real er12 = -0.2235530786388629525884427845e-01;

        static const real a21 = 5.26001519587677318785587544488e-2;

        static const real a31 = 1.97250569845378994544595329183e-2;
        static const real a32 = 5.91751709536136983633785987549e-2;

        static const real a41 = 2.95875854768068491816892993775e-2;
        static const real a43 = 8.87627564304205475450678981324e-2;

        static const real a51 = 2.41365134159266685502369798665e-1;
        static const real a53 = -8.84549479328286085344864962717e-1;
        static const real a54 = 9.24834003261792003115737966543e-1;

        static const real a61 = 3.7037037037037037037037037037e-2;
        static const real a64 = 1.70828608729473871279604482173e-1;
        static const real a65 = 1.25467687566822425016691814123e-1;

        static const real a71 = 3.7109375e-2;
        static const real a74 = 1.70252211019544039314978060272e-1;
        static const real a75 = 6.02165389804559606850219397283e-2;
        static const real a76 = -1.7578125e-2;

        static const real a81 = 3.70920001185047927108779319836e-2;
        static const real a84 = 1.70383925712239993810214054705e-1;
        static const real a85 = 1.07262030446373284651809199168e-1;
        static const real a86 = -1.53194377486244017527936158236e-2;
        static const real a87 = 8.27378916381402288758473766002e-3;

        static const real a91 = 6.24110958716075717114429577812e-1;
        static const real a94 = -3.36089262944694129406857109825e0;
        static const real a95 = -8.68219346841726006818189891453e-1;
        static const real a96 = 2.75920996994467083049415600797e1;
        static const real a97 = 2.01540675504778934086186788979e1;
        static const real a98 = -4.34898841810699588477366255144e1;

        static const real a101 = 4.77662536438264365890433908527e-1;
        static const real a104 = -2.48811461997166764192642586468e0;
        static const real a105 = -5.90290826836842996371446475743e-1;
        static const real a106 = 2.12300514481811942347288949897e1;
        static const real a107 = 1.52792336328824235832596922938e1;
        static const real a108 = -3.32882109689848629194453265587e1;
        static const real a109 = -2.03312017085086261358222928593e-2;

        static const real a111 = -9.3714243008598732571704021658e-1;
        static const real a114 = 5.18637242884406370830023853209e0;
        static const real a115 = 1.09143734899672957818500254654e0;
        static const real a116 = -8.14978701074692612513997267357e0;
        static const real a117 = -1.85200656599969598641566180701e1;
        static const real a118 = 2.27394870993505042818970056734e1;
        static const real a119 = 2.49360555267965238987089396762e0;
        static const real a1110 = -3.0467644718982195003823669022e0;

        static const real a121 = 2.27331014751653820792359768449e0;
        static const real a124 = -1.05344954667372501984066689879e1;
        static const real a125 = -2.00087205822486249909675718444e0;
        static const real a126 = -1.79589318631187989172765950534e1;
        static const real a127 = 2.79488845294199600508499808837e1;
        static const real a128 = -2.85899827713502369474065508674e0;
        static const real a129 = -8.87285693353062954433549289258e0;
        static const real a1210 = 1.23605671757943030647266201528e1;
        static const real a1211 = 6.43392746015763530355970484046e-1;

        if (dydt_in)
        {
            // copy dydt
            for (i = 0; i < n; i++)
                k1[i] = dydt_in[i];
        }
        else
        {
            // first correction
            for (i = 0; i < n; i++)
                yt[i] = y[i];
            ode->function(t, yt, k1);
        }

        // 2. corection
        for (i = 0; i < n; i++)
            yt[i] = y[i] + h * (a21 * k1[i]);
        ode->function(t + c2 * h, yt, k2);

        // 3. corection
        for (i = 0; i < n; i++)
            yt[i] = y[i] + h * (a31 * k1[i] + a32 * k2[i]);
        ode->function(t + c3 * h, yt, k3);

        // 4. corection
        for (i = 0; i < n; i++)
            yt[i] = y[i] + h * (a41 * k1[i] + a43 * k3[i]);
        ode->function(t + c4 * h, yt, k4);

        // 5. corection
        for (i = 0; i < n; i++)
            yt[i] = y[i] + h * (a51 * k1[i] + a53 * k3[i] + a54 * k4[i]);
        ode->function(t + c5 * h, yt, k5);

        // 6. corection
        for (i = 0; i < n; i++)
            yt[i] = y[i] + h * (a61 * k1[i] + a64 * k4[i] + a65 * k5[i]);
        ode->function(t + c6 * h, yt, k6);

        // 7. corection
        for (i = 0; i < n; i++)
            yt[i] = y[i] + h * (a71 * k1[i] + a74 * k4[i] + a75 * k5[i] + a76 * k6[i]);
        ode->function(t + c7 * h, yt, k7);

        // 8. corection
        for (i = 0; i < n; i++)
            yt[i] = y[i] + h * (a81 * k1[i] + a84 * k4[i] + a85 * k5[i] + a86 * k6[i] + a87 * k7[i]);
        ode->function(t + c8 * h, yt, k8);

        // 9. corection
        for (i = 0; i < n; i++)
            yt[i] = y[i] + h * (a91 * k1[i] + a94 * k4[i] + a95 * k5[i] + a96 * k6[i] + a97 * k7[i] + a98 * k8[i]);
        ode->function(t + c9 * h, yt, k9);

        // 10. corection
        for (i = 0; i < n; i++)
            yt[i] = y[i] + h * (a101 * k1[i] + a104 * k4[i] + a105 * k5[i] + a106 * k6[i] + a107 * k7[i] + a108 * k8[i] + a109 * k9[i]);
        ode->function(t + c10 * h, yt, k10);

        // 11. corection
        for (i = 0; i < n; i++)
            yt[i] = y[i] + h * (a111 * k1[i] + a114 * k4[i] + a115 * k5[i] + a116 * k6[i] + a117 * k7[i] + a118 * k8[i] + a119 * k9[i] + a1110 * k10[i]);
        ode->function(t + c11 * h, yt, k11);

        // 12. corection
        for (i = 0; i < n; i++)
            yt[i] = y[i] + h * (a121 * k1[i] + a124 * k4[i] + a125 * k5[i] + a126 * k6[i] + a127 * k7[i] + a128 * k8[i] + a129 * k9[i] + a1210 * k10[i] + a1211 * k11[i]);
        ode->function(t, yt, k12);

        // final value
        for (i = 0; i < n; i++)
        {
            k_help[i] = (b1 * k1[i] + b6 * k6[i] + b7 * k7[i] + b8 * k8[i] + b9 * k9[i] + b10 * k10[i] + b11 * k11[i] + b12 * k12[i]);
            y[i] = y[i] + h * k_help[i];
        }

        // calculate error
        real err3, err5;
        for (i = 0; i < n; i++)
        {
            err3 = (k_help[i] - bhh1 * k1[i] - bhh2 * k9[i] - bhh3 * k3[i]) * h;
            err5 = (er1 * k1[i] + er6 * k6[i] + er7 * k7[i] + er8 * k8[i] + er9 * k9[i] + er10 * k10[i] + er11 * k11[i] + er12 * k12[i]) * h;
            err[i] = err5*err5/sqrtl(0.01*err3*err3 + err5*err5);
        }

        if (dydt_out)
        {
            ode->function(t+h, y, dydt_out);
        }
    }

    int DoPr853::get_order() const
    {
        return 8;
    }

    int DoPr853::get_err_order() const
    {
        return 9;
    }
}

