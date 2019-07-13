#include <iostream>
#include "cavityFlow.h"

using namespace std;

int nx, ny , nt, nit, c, rho;
float dx, dy, nu, dt;


int main()
{
    cavityFlow flow;
    flow.setVariables();

    float **u, **v, **p, **b;
    u = flow.setInitialConditions(ny,nx);
    v = flow.setInitialConditions(ny, nx);
    p = flow.setInitialConditions(ny, nx);
    b = flow.setInitialConditions(ny ,nx);

    float** u_result = flow.cavity_flow(nt, u, v, dt, dx, dy, p, rho, nu);


    cout << u_result[35][8];

    return 0;
}
