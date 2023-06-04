#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define nNodos 1000
#define W 0.0085
#define gamma 2.0
#define dt 0.0005
#define tMax 0.20
#define U_SIZE 8 * nNodos

double* linspace(double start, double end, int num) {
    double step = (end - start) / (num - 1);
    double* array = malloc(num * sizeof(double));
    for (int i = 0; i < num; i++) {
        array[i] = start + i * step;
    }
    return array;
}

double pow2(double a) {
    return a * a;
}

double* getU(
    double* Bx, double* By, double* Bz, double* rho, double* p,
    double* Vx, double* Vy, double* Vz)
{
    double Mx, My, Mz, V, B, E;

    double* U = malloc(U_SIZE * sizeof(double));

    //parallel
    for (size_t i = 0; i < nNodos; i++)
    {
        Mx = rho[i] * Vx[i];
        My = rho[i] * Vy[i];
        Mz = rho[i] * Vz[i];

        V = sqrt(pow2(Vx[i]) + pow2(Vy[i]) + pow2(Vz[i]));
        B = sqrt(pow2(Bx[i]) + pow2(By[i]) + pow2(Bz[i]));
        E = p[i] / (gamma - 1) + 0.5 * rho[i] * (pow2(V)) + (pow2(B)) / 2;

        U[0 * nNodos + i] = Bx[i];
        U[1 * nNodos + i] = By[i];
        U[2 * nNodos + i] = Bz[i];
        U[3 * nNodos + i] = rho[i];
        U[4 * nNodos + i] = Mx;
        U[5 * nNodos + i] = My;
        U[6 * nNodos + i] = Mz;
        U[7 * nNodos + i] = E;
    }

    return U;
}

void calculateAllValues(
    int length,
    double* x,     
    double* rho_p,
    double* Vx_p,
    double* Vy_p,
    double* Vz_p,
    double* p_p,
    double* Bx_p,
    double* By_p,
    double* Bz_p ) {

    double rho1 = 1.0;
    double Vx1 = 0.0;
    double Vy1 = 0.0;
    double Vz1 = 0.0;
    double p1 = 1.0;
    double Bx1 = 0.75;
    double By1 = 1.0;
    double Bz1 = 0.0;

    double rho2 = 0.125;
    double Vx2 = 0.0;
    double Vy2 = 0.0;
    double Vz2 = 0.0;
    double p2 = 0.1;
    double Bx2 = 0.75;
    double By2 = -1.0;
    double Bz2 = 0.0;

    //paralell
    for (int i = 0; i < length; i++) {
        double tanhVal = tanh(x[i] / W);
            rho_p[i] = (rho2 + rho1) / 2 + ((rho2 - rho1) / 2) * tanhVal;
            Vx_p[i] = (Vx2 + Vx1) / 2 + ((Vx2 - Vx1) / 2) * tanhVal;
            Vy_p[i] = (Vy2 + Vy1) / 2 + ((Vy2 - Vy1) / 2) * tanhVal;
            Vz_p[i] = (Vz2 + Vz1) / 2 + ((Vz2 - Vz1) / 2) * tanhVal;
            p_p[i] = (p2 + p1) / 2 + ((p2 - p1) / 2) * tanhVal;
            Bx_p[i] = (Bx2 + Bx1) / 2 + ((Bx2 - Bx1) / 2) * tanhVal;
            By_p[i] = (By2 + By1) / 2 + ((By2 - By1) / 2) * tanhVal;
            Bz_p[i] = (Bz2 + Bz1) / 2 + ((Bz2 - Bz1) / 2) * tanhVal;
    }
}

int main() {
    double* x = linspace(-1.0, 1.0, nNodos);
    double dx = fabs(x[1] - x[0]);

    double* rho = malloc(nNodos * sizeof(double));
    double* Vx = malloc(nNodos * sizeof(double));
    double* Vy = malloc(nNodos * sizeof(double));
    double* Vz = malloc(nNodos * sizeof(double));
    double* p = malloc(nNodos * sizeof(double));
    double* Bx = malloc(nNodos * sizeof(double));
    double* By = malloc(nNodos * sizeof(double));
    double* Bz = malloc(nNodos * sizeof(double));

    calculateAllValues(
        nNodos,
        x,     
        rho,
        Vx,
        Vy,
        Vz,
        p,
        Bx,
        By,
        Bz
    );

    double* U = getU(Bx, By, Bz, rho, p, Vx, Vy, Vz);

    for (size_t i = 0; i < U_SIZE; i++)
    {
        printf("U[%ld] = %f\n", i, U[i]);
    }
    

    // Print the results
    // for (int i = 0; i < nNodos; i++) {
        // printf("rho[%d] = %f\n", i, rho[i]);
        // printf("Vx[%d] = %f\n", i, Vx[i]);
        // printf("Vy[%d] = %f\n", i, Vy[i]);
        // printf("Vz[%d] = %f\n", i, Vz[i]);
        // printf("p[%d] = %f\n", i, p[i]);
        // printf("Bx[%d] = %f\n", i, Bx[i]);
        // printf("By[%d] = %f\n", i, By[i]);
        // printf("Bz[%d] = %f\n", i, Bz[i]);
    // }

    free(x);
    free(rho);
    free(Vx);
    free(Vy);
    free(Vz);
    free(p);
    free(Bx);
    free(By);
    free(Bz);

    return 0;
}
