#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define nNodos 1000
#define W 0.0085
#define gamma 2.0
#define dt 0.0005
#define tMax 0.20

typedef enum Parameter {
    RHO, VX, VY, VZ,
    P, BX, BY, BZ
} Parameter_t;

double* linspace(double start, double end, int num) {
    double step = (end - start) / (num - 1);
    double* array = malloc(num * sizeof(double));
    for (int i = 0; i < num; i++) {
        array[i] = start + i * step;
    }
    return array;
}

void calculateValues(double* x, int length, double* result, Parameter_t param) {
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

    for (int i = 0; i < length; i++) {
        double tanhVal = tanh(x[i] / W);
        switch (param)
        {
        case RHO:
            result[i] = (rho2 + rho1) / 2 + ((rho2 - rho1) / 2) * tanhVal;
            break;
        case VX:
            result[i] = (Vx2 + Vx1) / 2 + ((Vx2 - Vx1) / 2) * tanhVal;
            break;
        case VY:
            result[i] = (Vy2 + Vy1) / 2 + ((Vy2 - Vy1) / 2) * tanhVal;
            break;
        case VZ:
            result[i] = (Vz2 + Vz1) / 2 + ((Vz2 - Vz1) / 2) * tanhVal;
            break;
        case P:
            result[i] = (p2 + p1) / 2 + ((p2 - p1) / 2) * tanhVal;
            break;
        case BX:
            result[i] = (Bx2 + Bx1) / 2 + ((Bx2 - Bx1) / 2) * tanhVal;
            break;
        case BY:
            result[i] = (By2 + By1) / 2 + ((By2 - By1) / 2) * tanhVal;
            break;
        case BZ:
            result[i] = (Bz2 + Bz1) / 2 + ((Bz2 - Bz1) / 2) * tanhVal;
            break;
        default:
            break;
        }
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

    calculateValues(x, nNodos, rho, RHO);
    calculateValues(x, nNodos, Vx, VX);
    calculateValues(x, nNodos, Vy, VY);
    calculateValues(x, nNodos, Vz, VZ);
    calculateValues(x, nNodos, p, P);
    calculateValues(x, nNodos, Bx, BX);
    calculateValues(x, nNodos, By, BY);
    calculateValues(x, nNodos, Bz, BZ);

    // Print the results
    for (int i = 0; i < nNodos; i++) {
        // printf("rho[%d] = %f\n", i, rho[i]);
        // printf("Vx[%d] = %f\n", i, Vx[i]);
        // printf("Vy[%d] = %f\n", i, Vy[i]);
        // printf("Vz[%d] = %f\n", i, Vz[i]);
        printf("p[%d] = %f\n", i, p[i]);
        // printf("Bx[%d] = %f\n", i, Bx[i]);
        // printf("By[%d] = %f\n", i, By[i]);
        // printf("Bz[%d] = %f\n", i, Bz[i]);
    }

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
