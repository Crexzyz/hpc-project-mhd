#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <float.h>
#include <unistd.h>
#include <string.h>
#include <stdbool.h>
#include "timer.h"
#include <omp.h>

#define NODES 1000
#define ROWS 8
#define W 0.0085
#define GAMMA 2.0
#define DT 0.0005
#define T_MAX 0.20
#define U_SIZE ROWS * NODES

void print_matrix(double* flat_matrix) {
    for(size_t row = 0; row < 8; ++row) {
        for(size_t col = 0; col < NODES; ++col) {
            printf("%.4f\t", flat_matrix[row * NODES + col]);
        }
        printf("\n");
    }
}

bool is_matrix_nan(double* flat_matrix) {
    for(size_t row = 0; row < 8; ++row) {
        for(size_t col = 0; col < NODES; ++col) {
            if(isnan(flat_matrix[row * NODES + col])) {
                return true;
            }
        }
    }

    return false;
}

double* linspace(double start, double end, int num) {
    double step = (end - start) / (num - 1);
    double* array = malloc(num * sizeof(double));
    for (int i = 0; i < num; i++) {
        array[i] = start + i * step;
    }
    return array;
}

void set_row_values(
    double* x, double* rho, double* Vx, double* Vy,
    double* Vz, double* p, double* Bx, double* By,
    double* Bz) {
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

    for (size_t i = 0; i < NODES; i++)
    {
        double tanh_val = tanh(x[i] / W);
        rho[i] = (rho2 + rho1) / 2 + ((rho2 - rho1) / 2) * tanh_val;
        Vx[i] = (Vx2 + Vx1) / 2 + ((Vx2 - Vx1) / 2) * tanh_val;
        Vy[i] = (Vy2 + Vy1) / 2 + ((Vy2 - Vy1) / 2) * tanh_val;
        Vz[i] = (Vz2 + Vz1) / 2 + ((Vz2 - Vz1) / 2) * tanh_val;
        p[i] = (p2 + p1) / 2 + ((p2 - p1) / 2) * tanh_val;
        Bx[i] = (Bx2 + Bx1) / 2 + ((Bx2 - Bx1) / 2) * tanh_val;
        By[i] = (By2 + By1) / 2 + ((By2 - By1) / 2) * tanh_val;
        Bz[i] = (Bz2 + Bz1) / 2 + ((Bz2 - Bz1) / 2) * tanh_val;
    }
}

double* getU(
    double* Bx, double* By, double* Bz,
    double* rho, double* p,
    double* Vx, double* Vy, double* Vz) {
    
    double* U = malloc(U_SIZE * sizeof(double));

    double* u_Bx = U + 0 * NODES;
    double* u_By = U + 1 * NODES;
    double* u_Bz = U + 2 * NODES;
    double* u_rho = U + 3 * NODES;
    double* Mx = U + 4 * NODES;
    double* My = U + 5 * NODES;
    double* Mz = U + 6 * NODES;
    double* E = U + 7 * NODES;

    for(size_t i = 0; i < NODES; ++ i) {
        u_Bx[i] = Bx[i];
        u_By[i] = By[i];
        u_Bz[i] = Bz[i];
        u_rho[i] = rho[i];
        Mx[i] = rho[i] * Vx[i];
        My[i] = rho[i] * Vy[i];
        Mz[i] = rho[i] * Vz[i];

        double v = sqrt(pow(Vx[i], 2) + pow(Vy[i], 2) + pow(Vz[i], 2));
        double b = sqrt(pow(Bx[i], 2) + pow(By[i], 2) + pow(Bz[i], 2));
        E[i] = p[i] / (GAMMA - 1.0) + 1.0/2 * rho[i] * pow(v, 2) + pow(b, 2) / 2;
    }

    return U;
}

double* getF(double* U) {
    double* F = malloc(U_SIZE * sizeof(double));
    memset(F, 0, U_SIZE * sizeof(double));

    double* Bx = U + 0 * NODES;
    double* By = U + 1 * NODES;
    double* Bz = U + 2 * NODES;
    double* rho = U + 3 * NODES;
    double* Mx = U + 4 * NODES;
    double* My = U + 5 * NODES;
    double* Mz = U + 6 * NODES;
    double* E = U + 7 * NODES;

    #pragma omp parallel for shared(F,Bx,By,Bz,rho,Mx,My,Mz,E)
    for (size_t i = 0; i < NODES; i++)
    {
        double Vx = Mx[i] / rho[i];
        double Vy = My[i] / rho[i];
        double Vz = Mz[i] / rho[i];

        F[0 * NODES + i] = 0;
        F[1 * NODES + i] = By[i]*Vx - Bx[i]*Vy;
        F[2 * NODES + i] = Bz[i]*Vx - Bx[i]*Vz;
        F[3 * NODES + i] = rho[i]*Vx;

        double b = sqrt(pow(Bx[i], 2) + pow(By[i], 2) + pow(Bz[i], 2));
        double v = sqrt(pow(Vx, 2) + pow(Vy, 2) + pow(Vz, 2));
        double p = (E[i] - 1.0/2 * rho[i] * pow(v, 2) - pow(b, 2)/2.0) * (GAMMA - 1.0);

        F[4 * NODES + i] = rho[i]*pow(Vx, 2) + (p + (pow(b, 2)/2.0)) - pow(Bx[i], 2);
        F[5 * NODES + i] = rho[i] * Vx * Vy - Bx[i] * By[i];
        F[6 * NODES + i] = rho[i] * Vx * Vz - Bx[i] * Bz[i];

        double F8a = ((1.0/2 * rho[i] * pow(v, 2) + (GAMMA*p / (GAMMA - 1.0)) + pow(b, 2)) * Vx);
        double F8b = (Bx[i]*(Bx[i]*Vx + By[i]*Vy + Bz[i]*Vz));

        F[7 * NODES + i] = F8a - F8b;
    }

    return F;
}

double* getD(double* U) {
    double* D = malloc(U_SIZE * sizeof(double));
    memset(D, 0, U_SIZE * sizeof(double));

    double* Bx = U + 0 * NODES;
    double* By = U + 1 * NODES;
    double* Bz = U + 2 * NODES;
    double* rho = U + 3 * NODES;
    double* Mx = U + 4 * NODES;
    double* My = U + 5 * NODES;
    double* Mz = U + 6 * NODES;
    double* E = U + 7 * NODES;

    #pragma omp parallel for shared(D,Bx,By,Bz,rho,Mx,My,Mz,E)
    for (size_t i = 0; i < NODES; i++) {
        double Vx = Mx[i] / rho[i];
        double Vy = My[i] / rho[i];
        double Vz = Mz[i] / rho[i];

        double b = sqrt(pow(Bx[i], 2) + pow(By[i], 2) + pow(Bz[i], 2));
        double v = sqrt(pow(Vx, 2) + pow(Vy, 2) + pow(Vz, 2));
        double p = (E[i] - 1.0/2 * rho[i] * pow(v, 2) - pow(b, 2)/2.0) * (GAMMA - 1.0);

        D[4 * NODES + i] = Vx;
        D[7 * NODES + i] = p / rho[i];
    }

    return D;
}

double* getLU(double* U, double* x) {
    double* LU = malloc(U_SIZE * sizeof(double));
    memset(LU, 0, U_SIZE * sizeof(double));

    double* F = getF(U);
    double* D = getD(U);

    double* dummy1 = malloc(U_SIZE * sizeof(double));
    memset(dummy1, 0, U_SIZE * sizeof(double));
    double* dummy2 = malloc(U_SIZE * sizeof(double));
    memset(dummy2, 0, U_SIZE * sizeof(double));

    double a1 = 1.0 / 60.0;
    double a2 = -3.0 / 20.0;
    double a3 = 3.0 / 4.0;
    double a4 = -3.0 / 4.0;
    double a5 = 3.0 / 20.0;
    double a6 = -1.0 / 60.0;


    //Lo intentamos :) y no funca 
    #pragma omp parallel shared(dummy1,F,a1,a2,a3,a4,a5,a6)
    {
        double arreglo[ROWS] = {0, 0, 0, 0, 0, 0, 0, 0};
        #pragma omp for 
        for (size_t i = 3; i < NODES - 3; i++) 
        {
            for (size_t row = 0; row < ROWS; row++)
            {
                arreglo[row] += a1*F[NODES*row + i+3];
                arreglo[row] += a2*F[NODES*row + i+2];
                arreglo[row] += a3*F[NODES*row + i+1];
                arreglo[row] += a4*F[NODES*row + i-1];
                arreglo[row] += a5*F[NODES*row + i-2];
                arreglo[row] += a6*F[NODES*row + i-3];
            }

            for (size_t row = 0; row < ROWS; row++) {
                dummy1[row * NODES + i] = arreglo[row];
            }
            memset(arreglo, 0, ROWS * sizeof(double));
        }
    }

    double b1 = 1.0 / 90.0;
    double b2 = -3.0 / 20.0;
    double b3 = 3.0 / 2.0;
    double b4 = -49.0 / 18.0;
    double b5 = 3.0 / 2.0;
    double b6 = -3.0 / 20.0;
    double b7 = 1.0 / 90.0;

    #pragma omp parallel shared(dummy2,D,b1,b2,b3,b4,b5,b6,b7)
    {
        double arreglo[ROWS] = {0, 0, 0, 0, 0, 0, 0, 0};
        #pragma omp for 
        for (size_t i = 3; i < NODES - 3; i++) {
            for (size_t row = 0; row < ROWS; row++)
            {
                arreglo[row] += b1*D[NODES*row + i+3];
                arreglo[row] += b2*D[NODES*row + i+2];
                arreglo[row] += b3*D[NODES*row + i+1];
                arreglo[row] += b4*D[NODES*row + i];
                arreglo[row] += b5*D[NODES*row + i-1];
                arreglo[row] += b6*D[NODES*row + i-2];
                arreglo[row] += b7*D[NODES*row + i-3];
            }

            for (size_t row = 0; row < ROWS; row++) {
                dummy2[row * NODES + i] = arreglo[row];
            }
            memset(arreglo, 0, ROWS * sizeof(double));
        }
    }

    double eta_v = 0.001;
    double eta_t = 0.001;
    double* rho = U + 3 * NODES;

    #pragma omp parallel for shared(dummy2,eta_v,eta_t,rho)
    for (size_t i = 0; i < NODES; i++) {
        dummy2[5 * NODES + i] = eta_v * rho[i] * dummy2[5 * NODES + i];
        dummy2[7 * NODES + i] = eta_t * rho[i] * dummy2[7 * NODES + i];
    }

    double dx = fabs(x[1] - x[0]);

    #pragma omp parallel for shared(LU,dummy1,dummy2,dx)
    for (size_t i = 0; i < U_SIZE; i++) {
        LU[i] = dummy1[i] / dx + dummy2[i] / pow(dx, 2);
    }

    free(dummy2);
    free(dummy1);
    free(D);
    free(F);
    return LU;
}

double* getNU(double* U, double* x){
    double* LU = getLU(U, x);

    double* U1fake = malloc(U_SIZE * sizeof(double));
    memset(U1fake, 0, U_SIZE * sizeof(double));

    #pragma omp parallel for shared(U1fake,U,LU)
    for (size_t i = 0; i < U_SIZE; i++)
    {
        U1fake[i] = U[i] + DT/2 * LU[i];
    }
    double* LU1 = getLU(U1fake, x);

    double* U2fake = malloc(U_SIZE * sizeof(double));
    memset(U2fake, 0, U_SIZE * sizeof(double));
    #pragma omp parallel for shared(U2fake,U,LU1)
    for (size_t i = 0; i < U_SIZE; i++)
    {
        U2fake[i] = U[i] + DT/2 * LU1[i];
    }
    double* LU2 = getLU(U2fake, x);

    double* U3fake = malloc(U_SIZE * sizeof(double));
    memset(U3fake, 0, U_SIZE * sizeof(double));

    #pragma omp parallel for shared(U3fake,U,LU2)
    for (size_t i = 0; i < U_SIZE; i++)
    {
        U3fake[i] = U[i] + DT * LU2[i];
    }
    double* LU3 = getLU(U3fake, x);

    double* NU = malloc(U_SIZE * sizeof(double));
    memset(NU, 0, U_SIZE * sizeof(double));

    #pragma omp parallel for shared(NU,U,LU,LU1,LU2,LU3)
    for (size_t i = 0; i < U_SIZE; i++)
    {
        NU[i] = U[i] + DT * (1.0/6*LU[i] + 1.0/3*LU1[i] + 1.0/3*LU2[i] + 1.0/6*LU3[i]);
    }

    free(U3fake);
    free(LU3);
    free(U2fake);
    free(LU2);
    free(U1fake);
    free(LU1);
    free(LU);
    return NU;
}

int main() {
    struct timespec tstart;
    cpu_timer_start(&tstart);

    double* x = linspace(-1.0, 1.0, NODES);

    double* rho = malloc(NODES * sizeof(double));
    double* Vx = malloc(NODES * sizeof(double));
    double* Vy = malloc(NODES * sizeof(double));
    double* Vz = malloc(NODES * sizeof(double));
    double* p = malloc(NODES * sizeof(double));
    double* Bx = malloc(NODES * sizeof(double));
    double* By = malloc(NODES * sizeof(double));
    double* Bz = malloc(NODES * sizeof(double));

    set_row_values(x, rho, Vx, Vy, Vz, p, Bx, By, Bz);
    double* U = getU(Bx, By, Bz, rho, p, Vx, Vy, Vz);

    double* U1 = getNU(U, x);
    double* U2 = getNU(U1, x);
    double* U3 = getNU(U2, x);

    double* UNMenos3 = U;
    double* UNMenos2 = U1;
    double* UNMenos1 = U2;
    double* Un = U3;
    
    double* LUnMenos3 = getLU(UNMenos3, x);
    double* LUnMenos2 = getLU(UNMenos2, x);
    double* LUnMenos1 = getLU(UNMenos1, x);

    double tiempo = 3*DT;
    int total_runs = 0;

    while(T_MAX > tiempo) {
        tiempo = tiempo + DT;
        double * LUn = getLU(Un, x);

        double* UIM = malloc(U_SIZE * sizeof(double));
        memset(UIM, 0, U_SIZE * sizeof(double));

        #pragma omp parallel for shared(UIM)
        for (size_t i = 0; i < U_SIZE; i++)
        {
            UIM[i] = Un[i] + DT*((55.0/24)*LUn[i] - (59.0/24)*LUnMenos1[i] + (37.0/24)*LUnMenos2[i] - (9.0/24)*LUnMenos3[i]);
        }
        double* LUIM = getLU(UIM, x);
        
        double* UNmasUno = malloc(U_SIZE * sizeof(double));
        memset(UNmasUno, 0, U_SIZE * sizeof(double));

        #pragma omp parallel for shared(UNmasUno)
        for (size_t i = 0; i < U_SIZE; i++)
        {
            UNmasUno[i] = Un[i] + DT*(9.0/24*LUIM[i] + 19.0/24*LUn[i] - 5.0/24*LUnMenos1[i] + 1.0/24*LUnMenos2[i]);
        }

        // double tol = 1e-15;
        double error = 1.0;
        int n = 0;

        while (error > 20*DBL_EPSILON)
        {
            double* LUNmasUno =  getLU(UNmasUno, x);

            double* UNmasUnoP = malloc(U_SIZE * sizeof(double));
            memset(UNmasUnoP, 0, U_SIZE * sizeof(double));

            #pragma omp parallel for shared(UNmasUnoP)
            for (size_t i = 0; i < U_SIZE; i++)
            {
                UNmasUnoP[i] = Un[i] + DT*((9.0/24)*LUNmasUno[i] + (19.0/24)*LUn[i] - (5.0/24)*LUnMenos1[i] + (1.0/24)*LUnMenos2[i]);
            }
            
            double max = -INFINITY;
            for (size_t i = 0; i < U_SIZE; i++) {
                double resta = UNmasUno[i] - UNmasUnoP[i];
                double local_error = fabs(resta);
                if(local_error > max) {
                    max = local_error;
                }
            }

            error = max;

            #pragma omp parallel for shared(UNmasUno)
            for (size_t i = 0; i < U_SIZE; i++) {
                UNmasUno[i] = UNmasUnoP[i];
            }

            n = n + 1;

            free(UNmasUnoP);
            free(LUNmasUno);
        }

        free(UNMenos3);
        UNMenos3 = UNMenos2;
        UNMenos2 = UNMenos1;
        UNMenos1 = Un;
        Un = UNmasUno;

        // print_matrix(Un);
        // exit(0);

        free(LUnMenos1);
        free(LUnMenos2);
        free(LUnMenos3);
        LUnMenos3 = getLU(UNMenos3,x);
        LUnMenos2 = getLU(UNMenos2,x);
        LUnMenos1 = getLU(UNMenos1,x);
        
        // plot(x,Un(4,:),'-b.')
        // pause(0.01)
        // sleep(0.01);
        total_runs += 1;

        free(LUIM);
        free(UIM);
        free(LUn);
    }

    // printf("%d\n", total_runs);
    // print_matrix(Un);

    free(UNMenos3);
    free(UNMenos2);
    free(UNMenos1);
    free(Un);
    free(LUnMenos1);
    free(LUnMenos2);
    free(LUnMenos3);
    free(Bz);
    free(By);
    free(Bx);
    free(p);
    free(Vz);
    free(Vy);
    free(Vx);
    free(rho);
    free(x);

    double elapsed = cpu_timer_stop(tstart);
    // printf("Total runs: %d\n", total_runs);
    printf("%f\n", elapsed);
    return 0;
}