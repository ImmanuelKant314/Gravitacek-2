#include "gravitacek2/geomotion/geomotion.hpp"

namespace gr2
{
    bool GeoMotion::necesarry_calculate(const real *y, real *y_save, const int& n)
    {
        if (!y_save)
        {
            y_save = new real[n];
            return true;
        }

        bool test = false;
        for (int i = 0; i < n; i++)
        {
            if (!test && y[i] != y_save[i])
                test = true;
            y_save[i] = y[i];
        }
        return test;
    }

    GeoMotion::GeoMotion(const int &dim, const int &n) : ODE(n)
    {
        // dimension
        this->dim = dim;

        // metric
        this->metric = new real*[dim];
        for (int i = 0; i < dim; i++)
        {
            metric[i] = new real[dim]{};
        }

        // christoffel symbols
        this->christoffel_symbols = new real**[dim];
        for (int i = 0; i < dim; i++)
        {
            this->christoffel_symbols[i] = new real*[dim];
            for (int j = 0; j < dim; j++)
            {
                this->christoffel_symbols[i][j] = new real[dim]{};
            }
        }

        // riemann tensor
        this->riemann_tensor = new real***[dim];
        for (int i = 0; i < dim; i++)
        {
            this->riemann_tensor[i] = new real**[dim];
            for (int j = 0; j < dim; j++)
            {
                this->riemann_tensor[i][j] = new real*[dim];
                for (int k = 0; k < dim; k++)
                {
                    this->riemann_tensor[i][j][k] = new real[dim]{};
                }
            }
        }

        // coordinates
        y_m = nullptr;
        y_c = nullptr;
        y_r = nullptr;
    };

    GeoMotion::~GeoMotion()
    {
        // metric
        for (int i = 0; i < dim; i++)
        {
            delete[] metric[i];
        }
        delete[] metric;

        // christoffel symbols
        for (int i = 0; i < dim; i++)
        {
            for (int j = 0; j < dim; j++)
                delete[] christoffel_symbols[i][j];
            delete[] christoffel_symbols[i];
        }
        delete[] christoffel_symbols;

        // riemann tensor
        for (int i = 0; i < dim; i++)
        {
            for (int j = 0; j < dim; j++)
            {
                for (int k = 0; k < dim; k++)
                    delete[] riemann_tensor[i][j][k];
                delete[] riemann_tensor[i][j];
            }
            delete[] riemann_tensor[i];
        }

        // coordinates
        delete[] y_c, y_m, y_r;
    };

    real **GeoMotion::get_metric() const
    {
        return metric;
    }

    real ***GeoMotion::get_christoffel_symbols() const
    {
        return christoffel_symbols;
    }

    real ****GeoMotion::get_riemann_tensor() const
    {
        return riemann_tensor;
    }

    void GeoMotion::function(const real &t, const real y[], real dydt[])
    {
        this->calculate_christoffel_symbols(y);
        
        // ========== Derivation of position ========== 
        for (int i = 0; i < dim; i++)
            dydt[i] = y[dim+i];

        // ========== Derivation of velocity ========== 
        for (int i = 0; i < dim; i++)
        {
            dydt[dim + i] = 0;
            for (int j = 0; j < dim; j++)
                for (int k = 0; k < dim; k++)
                    dydt[dim + i] += -christoffel_symbols[i][j][k]*y[dim+j]*y[dim+k];
        }
    }
}