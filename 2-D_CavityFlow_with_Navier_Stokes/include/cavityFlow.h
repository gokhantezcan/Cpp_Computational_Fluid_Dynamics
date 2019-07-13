#ifndef CAVITYFLOW_H
#define CAVITYFLOW_H


class cavityFlow
{
    private:
    protected:

    public:
        void setVariables();
        float** setInitialConditions(int ny, int nx);
        float** build_b(float** b, int rho, float dt, float** u, float** v, float dx, float dy);
        float** pressure_poisson(float** p, float dx, float dy, float** b);
        float**  cavity_flow(int nt, float** u, float** v, float dt, float dx, float dy, float** p, int rho, float nu);
};

#endif // CAVITYFLOW_H
