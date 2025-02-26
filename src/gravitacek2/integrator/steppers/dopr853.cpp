#include "gravitacek2/integrator/steppers.hpp"

#include <cmath>

namespace gr2
{
    DoPr853::DoPr853() : StepperBase(), k1(nullptr), k2(nullptr), k3(nullptr), k4(nullptr), k5(nullptr), k6(nullptr), k7(nullptr), k8(nullptr), k9(nullptr), k10(nullptr), k11(nullptr), k12(nullptr)
    {}

    DoPr853::~DoPr853()
    {
        delete[] k2, k3, k4, k5, k6, k7, k8, k9, k10, k11, k12;
        delete[] pc1, pc2, pc3, pc4, pc5, pc6, pc7, pc8;
        delete[] k_help;
    }

    void DoPr853::set_OdeSystem(std::shared_ptr<OdeSystem> ode)
    {
        int old_n = n;
        this->StepperBase::set_OdeSystem(ode);
        if(old_n != n)
        {
            delete[] k2, k3, k4, k5, k6, k7, k8, k9, k10, k11, k12;
            delete[] pc1, pc2, pc3, pc4, pc5, pc6, pc7, pc8;
            delete[] k_help;
            k1 = dydt_in;
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
            pc1 = new real[n];
            pc2 = new real[n];
            pc3 = new real[n];
            pc4 = new real[n];
            pc5 = new real[n];
            pc6 = new real[n];
            pc7 = new real[n];
            pc8 = new real[n];
        }
    }

    void DoPr853::reset()
    {
    }

    void DoPr853::step(const real &t, real y[], const real &h, const bool &dense, const real dydt_in[], real dydt_out[]) 
    {
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

         int i;
        
        // save time and step internaly
        this->t_in = t;
        this->h = h;

        // copy y to y_in, y_cur
        for (int i = 0; i < n; i++)
        {
            y_in[i] = y[i];
            y_cur[i] = y[i];
        }

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
                y_cur[i] = y_in[i];
            ode->function(t, y_cur, k1);
        }

        // 2. corection
        for (i = 0; i < n; i++)
            y_cur[i] = y[i] + h * (a21 * k1[i]);
        ode->function(t + c2 * h, y_cur, k2);

        // 3. corection
        for (i = 0; i < n; i++)
            y_cur[i] = y[i] + h * (a31 * k1[i] + a32 * k2[i]);
        ode->function(t + c3 * h, y_cur, k3);

        // 4. corection
        for (i = 0; i < n; i++)
            y_cur[i] = y[i] + h * (a41 * k1[i] + a43 * k3[i]);
        ode->function(t + c4 * h, y_cur, k4);

        // 5. corection
        for (i = 0; i < n; i++)
            y_cur[i] = y[i] + h * (a51 * k1[i] + a53 * k3[i] + a54 * k4[i]);
        ode->function(t + c5 * h, y_cur, k5);

        // 6. corection
        for (i = 0; i < n; i++)
            y_cur[i] = y[i] + h * (a61 * k1[i] + a64 * k4[i] + a65 * k5[i]);
        ode->function(t + c6 * h, y_cur, k6);

        // 7. corection
        for (i = 0; i < n; i++)
            y_cur[i] = y[i] + h * (a71 * k1[i] + a74 * k4[i] + a75 * k5[i] + a76 * k6[i]);
        ode->function(t + c7 * h, y_cur, k7);

        // 8. corection
        for (i = 0; i < n; i++)
            y_cur[i] = y[i] + h * (a81 * k1[i] + a84 * k4[i] + a85 * k5[i] + a86 * k6[i] + a87 * k7[i]);
        ode->function(t + c8 * h, y_cur, k8);

        // 9. corection
        for (i = 0; i < n; i++)
            y_cur[i] = y[i] + h * (a91 * k1[i] + a94 * k4[i] + a95 * k5[i] + a96 * k6[i] + a97 * k7[i] + a98 * k8[i]);
        ode->function(t + c9 * h, y_cur, k9);

        // 10. corection
        for (i = 0; i < n; i++)
            y_cur[i] = y[i] + h * (a101 * k1[i] + a104 * k4[i] + a105 * k5[i] + a106 * k6[i] + a107 * k7[i] + a108 * k8[i] + a109 * k9[i]);
        ode->function(t + c10 * h, y_cur, k10);

        // 11. corection
        for (i = 0; i < n; i++)
            y_cur[i] = y[i] + h * (a111 * k1[i] + a114 * k4[i] + a115 * k5[i] + a116 * k6[i] + a117 * k7[i] + a118 * k8[i] + a119 * k9[i] + a1110 * k10[i]);
        ode->function(t + c11 * h, y_cur, k11);

        // 12. corection
        for (i = 0; i < n; i++)
            y_cur[i] = y[i] + h * (a121 * k1[i] + a124 * k4[i] + a125 * k5[i] + a126 * k6[i] + a127 * k7[i] + a128 * k8[i] + a129 * k9[i] + a1210 * k10[i] + a1211 * k11[i]);
        ode->function(t, y_cur, k12);

        // final value
        for (i = 0; i < n; i++)
        {
            k_help[i] = (b1 * k1[i] + b6 * k6[i] + b7 * k7[i] + b8 * k8[i] + b9 * k9[i] + b10 * k10[i] + b11 * k11[i] + b12 * k12[i]);
            y_out[i] = y[i] + h * k_help[i];
            y[i] = y_out[i];
        }

        if (dydt_out)
        {
            ode->function(t+h, y, dydt_out);
            if (dense)
                for (i = 0; i < n; i++)
                    this->dydt_out[i] = dydt_out[i];
        }
        else if (dense)
        {
            ode->function(t+h, y, this->dydt_out);
        }
    }

    void DoPr853::step_err(const real &t, real y[], const real &h, real err[], const bool& dense, const real dydt_in[], real dydt_out[])
    {
        static const real c2 = 0.526001519587677318785587544488e-01;
        static const real c3 = 0.789002279381515978178381316732e-01;
        static const real c4 = 0.118350341907227396726757197510e+00;
        static const real c5 = 0.281649658092772603273242802490e+00;
        static const real c6 = 0.333333333333333333333333333333e+00;
        static const real c7 = 0.25e+00;
        static const real c8 = 0.30769230769230769230769230769+00;
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
        
        int i;
        
        // save time and step internaly
        this->t_in = t;
        this->h = h;

        // copy y to y_in, y_cur
        for (int i = 0; i < n; i++)
        {
            y_in[i] = y[i];
            y_cur[i] = y[i];
        }

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
                y_cur[i] = y_in[i];
            ode->function(t, y_cur, k1);
        }

        // 2. corection
        for (i = 0; i < n; i++)
            y_cur[i] = y[i] + h * (a21 * k1[i]);
        ode->function(t + c2 * h, y_cur, k2);

        // 3. corection
        for (i = 0; i < n; i++)
            y_cur[i] = y[i] + h * (a31 * k1[i] + a32 * k2[i]);
        ode->function(t + c3 * h, y_cur, k3);

        // 4. corection
        for (i = 0; i < n; i++)
            y_cur[i] = y[i] + h * (a41 * k1[i] + a43 * k3[i]);
        ode->function(t + c4 * h, y_cur, k4);

        // 5. corection
        for (i = 0; i < n; i++)
            y_cur[i] = y[i] + h * (a51 * k1[i] + a53 * k3[i] + a54 * k4[i]);
        ode->function(t + c5 * h, y_cur, k5);

        // 6. corection
        for (i = 0; i < n; i++)
            y_cur[i] = y[i] + h * (a61 * k1[i] + a64 * k4[i] + a65 * k5[i]);
        ode->function(t + c6 * h, y_cur, k6);

        // 7. corection
        for (i = 0; i < n; i++)
            y_cur[i] = y[i] + h * (a71 * k1[i] + a74 * k4[i] + a75 * k5[i] + a76 * k6[i]);
        ode->function(t + c7 * h, y_cur, k7);

        // 8. corection
        for (i = 0; i < n; i++)
            y_cur[i] = y[i] + h * (a81 * k1[i] + a84 * k4[i] + a85 * k5[i] + a86 * k6[i] + a87 * k7[i]);
        ode->function(t + c8 * h, y_cur, k8);

        // 9. corection
        for (i = 0; i < n; i++)
            y_cur[i] = y[i] + h * (a91 * k1[i] + a94 * k4[i] + a95 * k5[i] + a96 * k6[i] + a97 * k7[i] + a98 * k8[i]);
        ode->function(t + c9 * h, y_cur, k9);

        // 10. corection
        for (i = 0; i < n; i++)
            y_cur[i] = y[i] + h * (a101 * k1[i] + a104 * k4[i] + a105 * k5[i] + a106 * k6[i] + a107 * k7[i] + a108 * k8[i] + a109 * k9[i]);
        ode->function(t + c10 * h, y_cur, k10);

        // 11. corection
        for (i = 0; i < n; i++)
            y_cur[i] = y[i] + h * (a111 * k1[i] + a114 * k4[i] + a115 * k5[i] + a116 * k6[i] + a117 * k7[i] + a118 * k8[i] + a119 * k9[i] + a1110 * k10[i]);
        ode->function(t + c11 * h, y_cur, k11);

        // 12. corection
        for (i = 0; i < n; i++)
            y_cur[i] = y[i] + h * (a121 * k1[i] + a124 * k4[i] + a125 * k5[i] + a126 * k6[i] + a127 * k7[i] + a128 * k8[i] + a129 * k9[i] + a1210 * k10[i] + a1211 * k11[i]);
        ode->function(t, y_cur, k12);

        // final value
        for (i = 0; i < n; i++)
        {
            k_help[i] = (b1 * k1[i] + b6 * k6[i] + b7 * k7[i] + b8 * k8[i] + b9 * k9[i] + b10 * k10[i] + b11 * k11[i] + b12 * k12[i]);
            y_out[i] = y[i] + h * k_help[i];
            y[i] = y_out[i];
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
            if (dense)
                for (i = 0; i < n; i++)
                    this->dydt_out[i] = dydt_out[i];
        }
        else if (dense)
        {
            ode->function(t+h, y, this->dydt_out);
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

    void DoPr853::prepare_dense()
    {
        static const real c14 = 0.1e+00;
        static const real c15 = 0.2e+00;
        static const real c16 = 0.777777777777777777777777777778e+00;

        static const real a141 = 5.61675022830479523392909219681e-2;
        static const real a147 = 2.53500210216624811088794765333e-1;
        static const real a148 = -2.46239037470802489917441475441e-1;
        static const real a149 = -1.24191423263816360469010140626e-1;
        static const real a1410 = 1.5329179827876569731206322685e-1;
        static const real a1411 = 8.20105229563468988491666602057e-3;
        static const real a1412 = 7.56789766054569976138603589584e-3;
        static const real a1413 = -8.298e-3;
        static const real a151 = 3.18346481635021405060768473261e-2;
        static const real a156 = 2.83009096723667755288322961402e-2;
        static const real a157 = 5.35419883074385676223797384372e-2;
        static const real a158 = -5.49237485713909884646569340306e-2;
        static const real a1511 = -1.08347328697249322858509316994e-4;
        static const real a1512 = 3.82571090835658412954920192323e-4;
        static const real a1513 = -3.40465008687404560802977114492e-4;
        static const real a1514 = 1.41312443674632500278074618366e-1;
        static const real a161 = -4.28896301583791923408573538692e-1;
        static const real a166 = -4.69762141536116384314449447206e0;
        static const real a167 = 7.68342119606259904184240953878e0;
        static const real a168 = 4.06898981839711007970213554331e0;
        static const real a169 = 3.56727187455281109270669543021e-1;
        static const real a1613 = -1.39902416515901462129418009734e-3;
        static const real a1614 = 2.9475147891527723389556272149e0;
        static const real a1615 = -9.15095847217987001081870187138e0;

        static const real d41 = -0.84289382761090128651353491142e+01;
        static const real d46 = 0.56671495351937776962531783590e+00;
        static const real d47 = -0.30689499459498916912797304727e+01;
        static const real d48 = 0.23846676565120698287728149680e+01;
        static const real d49 = 0.21170345824450282767155149946e+01;
        static const real d410 = -0.87139158377797299206789907490e+00;
        static const real d411 = 0.22404374302607882758541771650e+01;
        static const real d412 = 0.63157877876946881815570249290e+00;
        static const real d413 = -0.88990336451333310820698117400e-01;
        static const real d414 = 0.18148505520854727256656404962e+02;
        static const real d415 = -0.91946323924783554000451984436e+01;
        static const real d416 = -0.44360363875948939664310572000e+01;

        static const real d51 = 0.10427508642579134603413151009e+02;
        static const real d56 = 0.24228349177525818288430175319e+03;
        static const real d57 = 0.16520045171727028198505394887e+03;
        static const real d58 = -0.37454675472269020279518312152e+03;
        static const real d59 = -0.22113666853125306036270938578e+02;
        static const real d510 = 0.77334326684722638389603898808e+01;
        static const real d511 = -0.30674084731089398182061213626e+02;
        static const real d512 = -0.93321305264302278729567221706e+01;
        static const real d513 = 0.15697238121770843886131091075e+02;
        static const real d514 = -0.31139403219565177677282850411e+02;
        static const real d515 = -0.93529243588444783865713862664e+01;
        static const real d516 = 0.35816841486394083752465898540e+02;

        static const real d61 = 0.19985053242002433820987653617e+02;
        static const real d66 = -0.38703730874935176555105901742e+03;
        static const real d67 = -0.18917813819516756882830838328e+03;
        static const real d68 = 0.52780815920542364900561016686e+03;
        static const real d69 =-0.11573902539959630126141871134e+02;
        static const real d610 = 0.68812326946963000169666922661e+01;
        static const real d611 = -0.10006050966910838403183860980e+01;
        static const real d612 = 0.77771377980534432092869265740e+00;
        static const real d613 = -0.27782057523535084065932004339e+01;
        static const real d614 = -0.60196695231264120758267380846e+02;
        static const real d615 = 0.84320405506677161018159903784e+02;
        static const real d616 = 0.11992291136182789328035130030e+02;

        static const real d71 = -0.25693933462703749003312586129e+02;
        static const real d76 = -0.15418974869023643374053993627e+03;
        static const real d77 = -0.23152937917604549567536039109e+03;
        static const real d78 = 0.35763911791061412378285349910e+03;
        static const real d79 = 0.93405324183624310003907691704e+02;
        static const real d710 = -0.37458323136451633156875139351e+02;
        static const real d711 = 0.10409964950896230045147246184e+03;
        static const real d712 = 0.29840293426660503123344363579e+02;
        static const real d713 = -0.43533456590011143754432175058e+02;
        static const real d714 = 0.96324553959188282948394950600e+02;
        static const real d715 = -0.39177261675615439165231486172e+02;
        static const real d716 = -0.14972683625798562581422125276e+03;

        int i;

        // simple coefficiens
        for (i = 0; i < n; i++)
        {
            pc1[i] = y_in[i];
            real ydiff = y_out[i]-y_in[i];
            pc2[i] = ydiff;
            real bspl = h*dydt_in[i] - ydiff;
            pc3[i] = bspl;
            pc4[i] = ydiff - h*dydt_out[i] - bspl;
            pc5[i] = d41*k1[i]+d46*k6[i]+d47*k7[i]+d48*k8[i]+d49*k9[i]+d410*k10[i]+d411*k11[i]+d412*k12[i];
            pc6[i] = d51*k1[i]+d56*k6[i]+d57*k7[i]+d58*k8[i]+d59*k9[i]+d510*k10[i]+d511*k11[i]+d512*k12[i];
            pc7[i] = d61*k1[i]+d66*k6[i]+d67*k7[i]+d68*k8[i]+d69*k9[i]+d610*k10[i]+d611*k11[i]+d612*k12[i];
            pc8[i] = d71*k1[i]+d76*k6[i]+d77*k7[i]+d78*k8[i]+d79*k9[i]+d710*k10[i]+d711*k11[i]+d712*k12[i];
        }

        // evaluation of function 1
        for (i = 0; i < n; i++)
            y_cur[i] = y_in[i]+h*(a141*k1[i]+a147*k7[i]+a148*k8[i]+a149*k9[i]+a1410*k10[i]+a1411*k11[i]+a1412*k12[i]+a1413*dydt_out[i]);
        this->ode->function(t_in+c14*h, y_cur, k10);

        // evaluation of function 2
        for (i = 0; i < n; i++)
            y_cur[i]=y_in[i]+h*(a151*k1[i]+a156*k6[i]+a157*k7[i]+a158*k8[i]+a1511*k11[i]+a1512*k12[i]+a1513*dydt_out[i]+a1514*k10[i]);
        this->ode->function(t_in+c15*h, y_cur, k2);

        // evaluation of function 3
        for (i = 0; i < n; i++)
            y_cur[i]=y_in[i]+h*(a161*k1[i]+a166*k6[i]+a167*k7[i]+a168*k8[i]+a169*k9[i]+a1613*dydt_out[i]+a1614*k10[i]+a1615*k2[i]);
        this->ode->function(t_in+c16*h, y_cur, k3);

        for (i = 0; i < n; i++)
        {
            pc5[i]=h*(pc5[i]+d413*dydt_out[i]+d414*k10[i]+d415*k2[i]+d416*k3[i]);
            pc6[i]=h*(pc6[i]+d513*dydt_out[i]+d514*k10[i]+d515*k2[i]+d516*k3[i]);
            pc7[i]=h*(pc7[i]+d613*dydt_out[i]+d614*k10[i]+d615*k2[i]+d616*k3[i]);
            pc8[i]=h*(pc8[i]+d713*dydt_out[i]+d714*k10[i]+d715*k2[i]+d716*k3[i]);
        }
    }

    real DoPr853::dense_out(const int &i, const real &t)
    {
        real s = (t-t_in)/h;
        real s1 = 1.0-s;
        return pc1[i]+s*(pc2[i]+s1*(pc3[i]+s*(pc4[i]+s1*(pc5[i]+s*(pc6[i]+s1*(pc7[i]+s*pc8[i]))))));
    }
}

