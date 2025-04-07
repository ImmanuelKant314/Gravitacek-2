#pragma once

#include "gravitacek2/geomotion/weyl.hpp"
#include "gravitacek2/geomotion/spacetimes.hpp"
#include "gravitacek2/integrator/integrator.hpp"
#include "gravitacek2/integrator/odesystems.hpp"
#include "gravitacek2/chaos/linearized_evolution.hpp"

#include <stdexcept>
#include <iostream>
#include <algorithm>
#include <fstream>
#include <array>
#include <cmath>
#include <chrono>
#include <thread>

class DataRecord : public gr2::Event
{
protected:
    int n;

public:
    std::vector<std::vector<gr2::real>> data;

    DataRecord(int n) : gr2::Event(gr2::EventType::data)
    {
        this->n = n;
    }

    virtual gr2::real value(const gr2::real &t, const gr2::real &dt, const gr2::real y[], const gr2::real dydt[]) override
    {
        return 0;
    }

    virtual void apply(gr2::StepperBase* stepper, gr2::real &t, gr2::real &dt, gr2::real y[], gr2::real dydt[]) override
    {
        std::vector<gr2::real> record;
        record.push_back(t);
        for (int i=0; i < n; i++)
            record.push_back(y[i]);
        this->data.push_back(record);
    }
};

template<class T>
class StopOnDisk: public gr2::Event
{
protected:
    std::shared_ptr<T> spt;
public:
    bool poincare;
    std::vector<std::array<gr2::real, 2>> data;
    gr2::real z;
    StopOnDisk(std::shared_ptr<T> spt, gr2::real z, bool poincare=false) : gr2::Event(gr2::EventType::modyfing), z(z), poincare(poincare), data(), spt(spt)
    {

    }
    virtual gr2::real value(const gr2::real &t, const gr2::real &dt, const gr2::real y[], const gr2::real dydt[]) override
    {
        int sign = y[gr2::Weyl::UZ]>0?1:-1;
        return (y[gr2::Weyl::Z]+sign*this->z);
    }

    virtual void apply(gr2::StepperBase* stepper, gr2::real &t, gr2::real &dt, gr2::real y[], gr2::real dydt[]) override
    {
        if (poincare)
            data.push_back({y[gr2::Weyl::RHO],y[gr2::Weyl::URHO]});
        y[gr2::Weyl::Z]*=-1;
        spt->function(t, y, dydt);
    }
};

class StopOnDiskTwoParticles : public gr2::Event
{
protected:
    std::shared_ptr<gr2::GeoMotion> spt;
    int sign;
    int n;
    bool poincare;
public:
    std::vector<std::array<gr2::real, 2>> data;
    gr2::real z;
    StopOnDiskTwoParticles(std::shared_ptr<gr2::GeoMotion> spt, gr2::real z, bool poincare=false) : gr2::Event(gr2::EventType::modyfing), z(z), spt(spt), poincare(poincare)
    {
        this->n = spt->get_n();
    }

    virtual gr2::real value(const gr2::real &t, const gr2::real &dt, const gr2::real y[], const gr2::real dydt[]) override
    {
        int sign = y[gr2::Weyl::UZ]>0?1:-1;
        return (y[gr2::Weyl::Z]+sign*this->z);
    }

    virtual void apply(gr2::StepperBase* stepper, gr2::real &dt, gr2::real &t, gr2::real y[], gr2::real dydt[]) override
    {
        if (poincare)
            data.push_back({y[gr2::Weyl::RHO],y[gr2::Weyl::URHO]});
        y[gr2::Weyl::Z]*=-1;
        (y+n)[gr2::Weyl::Z] = 2*y[gr2::Weyl::Z] + (y+n)[gr2::Weyl::Z];
        spt->function(t, y, dydt);
        spt->function(t, y+n, dydt+n);
    }
};

class RegularizeApproach : public gr2::Event
{
public:
    gr2::real z_min, z_app, alpha, beta;
    RegularizeApproach(gr2::real z_min, gr2::real z_app, gr2::real alpha, gr2::real beta):gr2::Event(gr2::data, false), z_min(z_min), z_app(z_app), alpha(alpha), beta(beta)
    {}

    virtual gr2::real value(const gr2::real &t, const gr2::real &dt, const gr2::real y[], const gr2::real dydt[]) override
    {
        gr2::real zz = std::abs(std::abs(y[gr2::Weyl::Z])-z_app);
        if (zz>z_min && std::abs(dt*y[gr2::Weyl::UZ]) > alpha*zz && y[gr2::Weyl::Z]*y[gr2::Weyl::UZ]<0) 
            return 0;
        else
            return 1;
    }

    virtual void apply(gr2::StepperBase* stepper, gr2::real &t, gr2::real &dt, gr2::real y[], gr2::real dydt[]) override
    {
        // std::cout << "Approach" << std::endl;
        // std::cout << dt << " " << y[3] << " " << y[gr2::Weyl::UZ] <<std::endl;
        gr2::real zz = std::abs(std::abs(y[gr2::Weyl::Z])-z_app);
        dt = beta*std::abs(y[gr2::Weyl::Z]/y[gr2::Weyl::UZ]);
        // std::cout << dt << " " << y[3] << std::endl;
        // std::cout << "Approach end" << std::endl;
    }
};

class StopBeforeBlackHole : public gr2::Event
{
public:
    gr2::real rho_min;
    bool activated;
    gr2::real t;

    StopBeforeBlackHole(gr2::real rho_min): gr2::Event(gr2::EventType::data, true), rho_min(rho_min), activated(false), t(0)
    {

    }

    virtual gr2::real value(const gr2::real &t, const gr2::real &dt, const gr2::real y[], const gr2::real dydt[]) override
    {
        return y[gr2::Weyl::RHO] < rho_min?0:1;
    }

    virtual void apply(gr2::StepperBase* stepper, gr2::real &t, gr2::real &dt, gr2::real y[], gr2::real dydt[]) override
    {
        activated=true;
        this->t = t;
    }
};

template<class T>
class StopTooHighErrorE : public gr2::Event
{
public:
    std::shared_ptr<T> spt;
    gr2::real E, E_, eps;
    bool activated;
    gr2::real t;

    StopTooHighErrorE(std::shared_ptr<T> spt, gr2::real E, gr2::real eps):gr2::Event(gr2::EventType::data, true), spt(spt), E(E), eps(eps), t(0), activated(false)
    {}

    virtual gr2::real value(const gr2::real &t, const gr2::real &dt, const gr2::real y[], const gr2::real dydt[]) override
    {
        spt->calculate_metric(y);
        E_ = - spt->get_metric()[gr2::Weyl::T][gr2::Weyl::T]*y[gr2::Weyl::UT];
        return abs(E_-E)/E < eps?1:0;
    }

    virtual void apply(gr2::StepperBase* stepper, gr2::real &t, gr2::real &dt, gr2::real y[], gr2::real dydt[]) override
    {
        activated=true;
        this->t = t;
    }
};

template<class T>
class StopTooHighErrorL : public gr2::Event
{
public:
    std::shared_ptr<T> spt;
    gr2::real L, L_, eps;
    bool activated;
    gr2::real t;

    StopTooHighErrorL(std::shared_ptr<T> spt, gr2::real L, gr2::real eps):gr2::Event(gr2::EventType::data, true), spt(spt), L(L), eps(eps), activated(false), t(0)
    {}

    virtual gr2::real value(const gr2::real &t, const gr2::real &dt, const gr2::real y[], const gr2::real dydt[]) override
    {
        spt->calculate_metric(y);
        L_ = spt->get_metric()[gr2::Weyl::PHI][gr2::Weyl::PHI]*y[gr2::Weyl::UPHI];
        return abs(L_-L)/L < eps?1:0;
    }

    virtual void apply(gr2::StepperBase* stepper, gr2::real &t, gr2::real &dt, gr2::real y[], gr2::real dydt[]) override
    {
        activated=true;
        this->t = t;
    }
};

class RenormalizationOfSecondParticleWeyl : public gr2::Event
{
protected:
    gr2::real target_norm;
    std::shared_ptr<gr2::GeoMotion> spt;
    int renorm_index;
public:
    gr2::real log_norm;

    RenormalizationOfSecondParticleWeyl(std::shared_ptr<gr2::GeoMotion> spt, gr2::real target_norm):gr2::Event(gr2::EventType::data, false), spt(spt), target_norm(target_norm), log_norm(0)
    {}

    virtual gr2::real value(const gr2::real &t, const gr2::real &dt, const gr2::real y[], const gr2::real dydt[]) override
    {
        return 0;
    }

    virtual void apply(gr2::StepperBase* stepper, gr2::real &t, gr2::real &dt, gr2::real y[], gr2::real dydt[]) override
    {
        // Calculate normalization 
        int n = spt->get_n();
        gr2::real norm2 = 0;
        gr2::real *y_ = y+n;
        for (int i = 0; i < 8; i++)
        {
            gr2::real dy = y[i]-y_[i];
            norm2 += dy*dy;
        }
        gr2::real norm = sqrtl(norm2);
        gr2::real factor = target_norm/norm;

        // Change y (pos, vels, lambda)
        for (int i = 0; i < n; i++)
            y_[i] = (y_[i]-y[i])*factor + y[i];

        // Change dydt ()
        spt->function(t, y_, dydt+n);

        // Save normalization
        this->log_norm += -log(factor);
    }
};

template <class T, int N>
class NumericalExpansions : public gr2::Event
{
public:
    gr2::real y_[N];
    gr2::real t_prev;
    gr2::real log_norm_prev;
    gr2::real total_norm_prev;

    gr2::real t_last_step;

    std::shared_ptr<T> spt;
    gr2::real rho_min, rho_max, n_rho;
    gr2::real delta_rho;
    gr2::real z_min, z_max, n_z;
    gr2::real delta_z;
    gr2::real *log_norm;
    gr2::real dt;
    gr2::real** data;
    gr2::real** time_spend_in_area;

    bool test;

    NumericalExpansions(std::shared_ptr<T> spt, gr2::real rho_min, gr2::real rho_max, int n_rho, gr2::real z_min, gr2::real z_max, int n_z, gr2::real *log_norm):gr2::Event(gr2::EventType::data, false), spt(spt), dt(dt), rho_min(rho_min), rho_max(rho_max), n_rho(n_rho), z_min(z_min), z_max(z_max), n_z(n_z), log_norm(log_norm), t_prev(0), test(false)
    {
        delta_rho = (rho_max-rho_min)/n_rho;
        delta_z = (z_max-z_min)/n_z;
        data = new gr2::real*[n_rho];
        time_spend_in_area = new gr2::real*[n_rho];
        for(int i = 0; i < n_rho; i++)
        {
            data[i] = new gr2::real[n_z]{};
            time_spend_in_area[i] = new gr2::real[n_z]{};
        }

        // time of previsou step
        t_prev = 0;
        t_last_step = 0;
    }

    ~NumericalExpansions()
    {
        for (int i = 0; i < n_rho; i++)
        {
            delete[] data[i];
            delete[] time_spend_in_area[i];
        }
        delete[] data;
        delete[] time_spend_in_area;
    }

    virtual gr2::real value(const gr2::real &t, const gr2::real &dt, const gr2::real y[], const gr2::real dydt[]) override
    {
        return 0;
    }

    virtual void apply(gr2::StepperBase* stepper, gr2::real &t, gr2::real &dt_, gr2::real y[], gr2::real dydt[]) override
    {
        // std::cout << "========== Apply ==========" << std::endl;
        // TODO: for-cycle (check if any line was crossed)
        // iterate from prev. time to current time
        // check both z and rho
        // than do linear interpolation of values
        // this has to be done before renormalization
        // indices
        int n = spt->get_n();
        int dim = spt->get_dim();

        int rho_index = gr2::Weyl::RHO;
        int z_index = gr2::Weyl::Z;
        int urho_index = gr2::Weyl::URHO;
        int uz_index = gr2::Weyl::UZ;

        // prepare iteration
        gr2::real urho = std::max(abs(stepper->dense_out(urho_index, t_last_step)), abs(stepper->dense_out(rho_index, t)));
        gr2::real uz = std::max(abs(stepper->dense_out(uz_index, t_last_step)), abs(stepper->dense_out(uz_index, t)));
        gr2::real dt = std::min(delta_rho, delta_z)/sqrtl(urho*urho + uz*uz)/10; // estimate time step for rho direction
        int iters = int((t-t_last_step)/dt) + 2;

        // iterate
        gr2::real rho_iter_old = stepper->dense_out(rho_index, t_last_step);
        gr2::real z_iter_old = stepper->dense_out(z_index, t_last_step);
        gr2::real rho_iter_new, z_iter_new;
        int i_iter_old = int((rho_iter_old-rho_min)/delta_rho);
        int j_iter_old = int((z_iter_old-z_min)/delta_z);
        int i_iter_new, j_iter_new;

        gr2::real t_step = (t-t_last_step)/(iters-1);

        // std::cout << "Brum" << std::endl;
        for (int i = 1; i < iters; i++)
        {
            gr2::real t_now = t_last_step+t_step*i;

            // new position
            rho_iter_new = stepper->dense_out(rho_index, t_now);
            z_iter_new = stepper->dense_out(z_index, t_now);

            // new indices
            i_iter_new = int((rho_iter_new-rho_min)/delta_rho);
            j_iter_new = int((z_iter_new-z_min)/delta_z);

            int counter_of_events = 0;
            gr2::real t_event; 
            gr2::real t_event_rho = -1, t_event_z = -1;

            // calculate rho event
            if (i_iter_new - i_iter_old != 0)
            {
                gr2::real rho_target = std::max(i_iter_new, i_iter_old)*delta_rho+rho_min;
                // calculate time
                t_event_rho = t_now - t_step*(rho_target-rho_iter_new)/(rho_iter_old-rho_iter_new);
                counter_of_events++;
            }

            // calculate z event
            if (j_iter_new - j_iter_old != 0)
            {
                gr2::real z_target = std::max(j_iter_new, j_iter_old)*delta_z+z_min;
                // calculate time
                t_event_rho = t_now - t_step*(z_target-z_iter_new)/(z_iter_old-z_iter_new);
                counter_of_events++;
            }

            // take the last
            if (counter_of_events != 0)
            {
                // calculate time
                t_event = std::max(t_event_rho, t_event_z);

                // get position
                for (int j = 0; j < N; j++)
                    y_[j] = stepper->dense_out(j, t_event);
                spt->calculate_metric(y_);
                spt->calculate_christoffel_symbols(y_);

                // paralel transport 
                for (int j = 0; j < spt->get_dim(); j++)
                    for (int k = 0; k < spt->get_dim(); k++)
                        for (int l = 0; l < spt->get_dim(); l++)
                            y_[n + dim + j] += spt->get_christoffel_symbols()[j][k][l]*y_[dim+k]*(y_[n+l] - y_[l]);

                // projection
                gr2::real u_up_indices[4]{};
                gr2::real u_down_indices[4]{};

                for (int j = 0; j < 4; j++)
                    u_up_indices[j] = y_[4+j];

                for (int j = 0; j < 4; j++)
                    for (int k = 0; k < 4; k++)
                        u_down_indices[j] += spt->get_metric()[j][k]*u_up_indices[k];

                // calculate norm of separation in space 
                gr2::real norm_of_sep2 = 0;
                gr2::real dyj, dyk, dvj, dvk;
                for (int j = 0; j < 4; j++)
                    for (int k = 0; k < 4; k++)
                    {
                        dyj = y_[n+j]-y_[j];
                        dyk = y_[n+k]-y_[k];
                        dvj = y_[n+dim+j]-y_[dim+j];
                        dvk = y_[n+dim+k]-y_[dim+k];
                        norm_of_sep2 += (spt->get_metric()[j][k]+u_down_indices[j]*u_down_indices[k])*dyj*dyk;
                        norm_of_sep2 += (spt->get_metric()[j][k]+u_down_indices[j]*u_down_indices[k])*dvj*dvk;
                    }
                gr2::real log_norm_of_sep = 0.5*log(norm_of_sep2);

                // save result to array
                // std::cout << (log_norm_of_sep + *(this->log_norm) - log_norm_prev) << std::endl;
                if (test)
                {
                    this->data[i_iter_old][j_iter_old] += log_norm_of_sep + *(this->log_norm) - log_norm_prev;
                    this->time_spend_in_area[i_iter_old][j_iter_old] += t_event - t_prev;
                }
                test = true;
                t_prev = t_event;
                log_norm_prev = log_norm_of_sep + *(this->log_norm);
            }

            // new to old
            rho_iter_old = rho_iter_new;
            z_iter_old = z_iter_new;
            i_iter_old = i_iter_new;
            j_iter_old = j_iter_new;
        }

        // save values for nex step
        t_last_step = t;
    }
};

class ConstantStepDataMonitoring : public gr2::Event
{
public:
    gr2::real t;
    gr2::real h;
    std::vector<gr2::real> times;
    std::vector<std::array<gr2::real, 9>> data;
    ConstantStepDataMonitoring(gr2::real t_init, gr2::real h) : gr2::Event(gr2::EventType::data), times(), data()
    {
        t = t_init;
        this->h = h;
    }
    virtual gr2::real value(const gr2::real &t, const gr2::real &dt, const gr2::real y[], const gr2::real dydt[]) override
    {   
        return 0;
    }
    virtual void apply(gr2::StepperBase* stepper, gr2::real &t, gr2::real &dt, gr2::real y[], gr2::real dydt[]) override
    {
        while (this->t<t)
        {
            times.push_back(this->t);
            std::array<gr2::real, 9> y_val;
            for (int i = 0; i<9; i++)
            {
                y_val[i] = stepper->dense_out(i, this->t);
            }
            data.push_back(y_val);
            this->t += h;
        }
    };
};