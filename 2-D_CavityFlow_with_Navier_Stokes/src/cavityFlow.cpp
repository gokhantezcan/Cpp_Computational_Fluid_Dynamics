#include "cavityFlow.h"
#include <iostream>

using namespace std;

extern int nx, ny, nt, nit, c, rho;
extern float dx, dy, nu, dt;

void cavityFlow::setVariables()
{
    nx = 41, ny = 41, nt = 500, nit = 50, c = 1, rho = 1;
    dx = 0.05, dy = 0.05, nu = 0.1, dt = 0.001;

}

float** cavityFlow::setInitialConditions(int ny, int nx)
{
    float **a;
    a = new float *[ny];
    for(int i = 0; i < ny; i++)
        a[i] = new float[nx];

    for(int j = 0; j < ny; j++)
    {
        for(int i = 0; i < nx; i++)
        {
            a[j][i] = 0.0;
        }
    }
    return a;
}

float** cavityFlow::build_b(float** b2, int rho, float dt, float** u, float** v, float dx, float dy)
{
    for(int j = 1; j < ny-1; j++)
    {
        for(int i = 1; i < nx-1; i++)
        {
            b2[j][i] = (((1000) * (((u[j][i+1] - u[j][i-1]) / (0.1)) + ((v[j+1][i] - v[j-1][i]) / (0.1)))) -
                    (((u[j][i+1] - u[j][i-1]) / (0.1)) * ((u[j][i+1] - u[j][i-1]) / (0.1))) -
                    (2 * ((u[j+1][i] - u[j-1][i]) / (0.1)) * ((v[j][i+1] - v[j][i-1]) / (0.1))) -
                    (((v[j+1][i] - v[j-1][i]) / (0.1)) * ((v[j+1][i] - v[j-1][i]) / (0.1))));
        }
    }

    return b2;
}

float** cavityFlow::pressure_poisson(float** p, float dx, float dy, float** b3)
{
    float** pn;
    pn = new float* [ny];
    for(int i = 0; i < ny; i++)
        pn[i] = new float[nx];

    for(int j = 0; j < ny; j++)
    {
        for(int i = 0; i < nx; i++)
        {
            pn[j][i] = p[j][i];
        }
    }


    for(int q = 0; q < nit; q++)
    {
        for(int j = 0; j < ny; j++)
        {
            for(int i = 0; i < nx; i++)
            {
                pn[j][i] = p[j][i];
            }
        }

        for(int j = 1; j < ny-1; j++)
        {
            for(int i = 1; i < nx-1; i++)
            {
                p[j][i] = (((((pn[j][i+1] + pn[j][i-1]) * (0.0025)) + ((pn[j+1][i] + pn[j-1][i]) * (0.0025))) / (0.01)) -
                        ((0.000625) * b3[j][i]));
            }
        }

        for(int j = 0; j < ny; j++)
        {
            p[j][40] = p[j][39];
            p[j][0] = p[j][1];
        }

        for(int i = 0; i < nx; i++)
        {
            p[0][i] = p[1][i];
            p[40][i] = 0.0;

        }

    }

    return p;
}

float**  cavityFlow::cavity_flow(int nt, float** u, float** v, float dt, float dx, float dy, float** p, int rho, float nu)
{
    float** un;
    float** vn;
    un = new float* [ny];
    vn = new float* [ny];
    for(int i = 0; i < ny; i++)
    {
        un[i] = new float[nx];
        vn[i] = new float[nx];
    }

    for(int j = 0; j < ny; j++)
    {
        for(int i = 0; i < nx; i++)
        {
            un[j][i] = u[j][i];
            vn[j][i] = v[j][i];
        }
    }

    float** b2;
    b2 = new float* [ny];
    for(int i = 0; i < ny; i++)
        b2[i] = new float[nx];
    for(int j = 0; j < ny; j++)
    {
        for(int i = 0; i < nx; i++)
        {
            b2[j][i] = 0.0;
        }
    }


    for(int n = 0; n < nit; n++)
    {
        for(int j = 0; j < ny; j++)
        {
            for(int i = 0; i < nx; i++)
            {
                un[j][i] = u[j][i];
                vn[j][i] = v[j][i];
            }
        }

        b2 = build_b(b2, rho, dt, u, v, dx, dy);
        p = pressure_poisson(p, dx, dy, b2);

        for(int j = 1; j < ny-1; j++)
        {
            for(int i = 1; i < nx-1; i++)
            {
                u[j][i] = (un[j][i] - ((un[j][i] * (0.02)) * (un[j][i] - un[j][i-1])) - ((vn[j][i] * (0.02)) * (un[j][i] - un[j-1][i])) -
                        ((0.01) * (p[j][i+1] - p[j][i-1])) + ((0.1) * ((0.4) * (un[j][i+1] - 2 * un[j][i] + un[j][i-1]) +
                        (0.4) * (un[j+1][i] - 2 * un[j][i] + un[j-1][i]))));



                v[j][i] = vn[j][i] - (un[j][i] * (0.02) * (vn[j][i] - vn[j][i-1])) - (vn[j][i] * (0.02) * (vn[j][i] - vn[j-1][i])) -
                        ((0.01) * (p[j+1][i] - p[j-1][i])) +  ((nu) * ((0.4) * (vn[j][i+1] - 2 * vn[j][i] + vn[j][i-1]) +
                        (0.4) * (vn[j+1][i] - 2 * vn[j][i] + vn[j-1][i])));

            }
        }
//
//
//
        for(int i = 0; i < nx; i++)
        {
            u[0][i] = 0.0;
            u[40][i] = 1.0;
            v[0][i] = 0.0;
            v[40][i] = 0.0;
        }

        for(int j = 0; j < ny; j++)
        {
            u[j][0] = 0.0;
            u[j][40] = 0.0;
            v[j][0] = 0.0;
            v[j][40] = 0.0;
        }
//
    }


    return u;
}

