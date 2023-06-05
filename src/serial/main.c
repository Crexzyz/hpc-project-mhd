#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <float.h>
#include <unistd.h>

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

double* getF(double* U)
{
    double* F = malloc(U_SIZE * sizeof(double));

    double* Bx = U + (0 * nNodos);
    double* By = U + (1 * nNodos);
    double* Bz = U + (2 * nNodos);
    double* rho = U + (3 * nNodos);
    double* Vx = U + (4 * nNodos);
    double* Vy = U + (5 * nNodos);
    double* Vz = U + (6 * nNodos);
    double* E = U + (7 * nNodos);

    double* B = malloc(nNodos * sizeof(double));
    double* V = malloc(nNodos * sizeof(double));
    double* p = malloc(nNodos * sizeof(double));

    for (size_t i = 0; i < nNodos; i++){
        Vx[i] = Vx[i] / rho[i];
        Vy[i] = Vy[i] / rho[i];
        Vz[i] = Vz[i] / rho[i];
    }

    //parallel
    for (size_t i = 0; i < nNodos; i++)
    {
        // gamma = 2.0;

        B[i] = sqrt(pow2(Bx[i]) + pow2(By[i]) + pow2(Bz[i]));
        V[i] = sqrt(pow2(Vx[i]) + pow2(Vy[i]) + pow2(Vz[i]));
        p[i] = (E[i] - 0.5 * rho[i] *  pow2(V[i]) - pow2(B[i]) / 2) * (gamma - 1);

        F[0 * nNodos + i] = 0;
        F[1 * nNodos + i] = By[i] * Vx[i] - Bx[i] * Vy[i];
        F[2 * nNodos + i] = Bz[i] * Vx[i] - Bx[i] * Vz[i];
        F[3 * nNodos + i] = rho[i] * Vx[i];
        F[4 * nNodos + i] = rho[i] * Vx[i] * Vx[i] + (p[i] + (pow2(B[i])) / 2) - pow2(Bx[i]);
        F[5 * nNodos + i] = rho[i] * Vx[i] * Vy[i] - Bx[i] * By[i];
        F[6 * nNodos + i] = rho[i] * Vx[i] * Vz[i] - Bx[i] * Bz[i];
        double F8a = ((0.5 * rho[i] * pow2(V[i]) + (gamma * p[i] / (gamma - 1)) + pow2(B[i])) * Vx[i]);
        double F8b = (Bx[i] * (Bx[i] * Vx[i] + By[i] * Vy[i] + Bz[i] * Vz[i]));
        F[7 * nNodos + i] = F8a - F8b;
    }

    free(p);
    free(V);
    free(B);
    return F;
}

double* getD(double* U)
{
    double* D = malloc(U_SIZE * sizeof(double));
    for (size_t i = 0; i < U_SIZE; i++)
    {
        D[i] = 0;
    }

    double* Bx = U + (0 * nNodos);
    double* By = U + (1 * nNodos);
    double* Bz = U + (2 * nNodos);
    double* rho = U + (3 * nNodos);
    double* Vx = U + (4 * nNodos);
    double* Vy = U + (5 * nNodos);
    double* Vz = U + (6 * nNodos);
    double* E = U + (7 * nNodos);

    double* B = malloc(nNodos * sizeof(double));
    double* V = malloc(nNodos * sizeof(double));
    double* p = malloc(nNodos * sizeof(double));

    for (size_t i = 0; i < nNodos; i++){
        Vx[i] = Vx[i] / rho[i];
        Vy[i] = Vy[i] / rho[i];
        Vz[i] = Vz[i] / rho[i];
    }

    for (size_t i = 0; i < nNodos; i++)
    {
        B[i] = sqrt(pow2(Bx[i]) + pow2(By[i]) + pow2(Bz[i]));
        V[i] = sqrt(pow2(Vx[i]) + pow2(Vy[i]) + pow2(Vz[i]));
        p[i] = (E[i] - 0.5 * rho[i] *  pow2(V[i]) - pow2(B[i]) / 2) * (gamma - 1);
        
        D[4 * nNodos + i] = Vx[i];
        D[7 * nNodos + i] = p[i] / rho[i];
    }

    free(p);
    free(V);
    free(B);
    return D;
}

void add_column(double* src_col, double* dst_mtx, size_t dst_col) {
    for(size_t row = 0; row < 8; ++row) {
        // printf("dst (%f) += %f\n", dst_mtx[row * nNodos + dst_col], src_col[row]);
        dst_mtx[row * nNodos + dst_col] += src_col[row];
    }
    // printf("End\n\n");
}

double* multiply_column(double value, size_t col, double *mtx) {
    double *a_result = malloc(8 * sizeof(double));

    for(size_t row = 0; row < 8; ++row) {
        a_result[row] = value * mtx[row * nNodos + col];
    }

   return a_result;
}

double* getLU(double* U, double *x) {
    double* rho = U + (3 * nNodos);
    double length = nNodos;
    double etaV = 0.001;
    double etaT = 0.001;
    double a1 = 1.0 / 60;
    double a2 = -3.0 / 20;
    double a3 = 3.0 / 4;
    double a4 = -3.0 / 4;
    double a5 = 3.0 / 20;
    double a6 = -1.0 / 60;
    
    double b1 = 1.0 / 90;
    double b2 = -3.0 / 20;
    double b3 = 3.0 / 2;
    double b4 = -49.0 / 18;
    double b5 = 3.0 / 2;
    double b6 = -3.0 / 20;
    double b7 = 1.0 / 90;
    
    // Allocate memory for arrays
    double* F = getF(U);
    double* D = getD(U);
    double* dummy1 = (double *)malloc(8 * length * sizeof(double));
    for (size_t i = 0; i < 8 * length; i++)
    {
        dummy1[i] = 0;
    }

    double* dummy2 = (double *)malloc(8 * length * sizeof(double));
    for (size_t i = 0; i < 8 * length; i++)
    {
        dummy2[i] = 0;
    }

    // Calculate dummy1 array
    //parallel
    for (size_t i = 3; i < length - 3; i++)
    {
        double* a1_val = multiply_column(a1, i+3, F);
        double* a2_val = multiply_column(a2, i+2, F);
        double* a3_val = multiply_column(a3, i+1, F);
        double* a4_val = multiply_column(a4, i-1, F);
        double* a5_val = multiply_column(a5, i-2, F);
        double* a6_val = multiply_column(a6, i-3, F);
        
        double arreglo[8] = {0, 0, 0, 0, 0, 0, 0, 0};
        for (size_t arr_i = 0; arr_i < 8; arr_i++){
            arreglo[arr_i] = a1_val[arr_i] + a2_val[arr_i] + a3_val[arr_i]+ a4_val[arr_i] + a5_val[arr_i] + a6_val[arr_i];
        }

        for(size_t row = 0; row < 8; ++row) {
            dummy1[row * nNodos + i] = arreglo[row];
        }

        free(a1_val);
        free(a2_val);
        free(a3_val);
        free(a4_val);
        free(a5_val);
        free(a6_val);
    }

    // Calculate dummy2 array
    for (size_t i = 3; i < length - 3; i++)
    {
        double* b1_val = multiply_column(b1, i + 3, D);
        double* b2_val = multiply_column(b2, i + 2, D);
        double* b3_val = multiply_column(b3, i + 1, D);
        double* b4_val = multiply_column(b4, i, D);
        double* b5_val = multiply_column(b5, i - 1, D);
        double* b6_val = multiply_column(b6, i - 2, D);
        double* b7_val = multiply_column(b7, i - 3, D);
        
        double arreglo[8] = {0, 0, 0, 0, 0, 0, 0, 0};

        for (size_t arr_i = 0; arr_i < 8; arr_i++){
            arreglo[arr_i] = b1_val[arr_i] + b2_val[arr_i] + b3_val[arr_i]+ b4_val[arr_i] + b5_val[arr_i] + b6_val[arr_i] + b7_val[arr_i];
        }   

        for(size_t row = 0; row < 8; ++row) {
            dummy2[row * nNodos + i] = arreglo[row];
        }

        free(b1_val);
        free(b2_val);
        free(b3_val);
        free(b4_val);
        free(b5_val);
        free(b6_val);
        free(b7_val);
    }

    // Apply additional calculations to dummy2 array
    for (size_t i = 0; i < length; i++)
    {
        dummy2[nNodos * 4 + i] = etaV * rho[i] * dummy2[nNodos * 4 + i];
        dummy2[nNodos * 7 + i] = etaT * rho[i] * dummy2[nNodos * 7 + i];
    }


    double* LU = malloc(U_SIZE * sizeof(double));
    for (size_t i = 0; i < U_SIZE; i++)
    {
        LU[i] = 0;
    }

    double dx = fabs(x[1] - x[0]);
    // Calculate LU array
    for (size_t i = 0; i < U_SIZE; i++)
    {
        LU[i] = (dummy1[i] / dx) + (dummy2[i] / pow2(dx));
        // printf("LU[%ld] = %.10e\n", i, LU[i]);
    }

    // Free allocated memory
    free(F);
    free(D);
    free(dummy1);
    free(dummy2);

    return LU;
}

double* getNU(double* U, double* x){
    double* LU = getLU(U, x);
    
    double* U1fake = malloc(U_SIZE * sizeof(double));
    for (size_t i = 0; i < U_SIZE; i++)
    {
        U1fake[i] = U[i] + dt/2*LU[i];
    }
    double* LU1 = getLU(U1fake, x);

    double* U2fake = malloc(U_SIZE * sizeof(double));
    for (size_t i = 0; i < U_SIZE; i++)
    {
        U2fake[i] = U[i] + dt/2 * LU1[i];
    }
    double* LU2 = getLU(U2fake, x);  

    double* U3fake = malloc(U_SIZE * sizeof(double));
    for (size_t i = 0; i < U_SIZE; i++)
    {
        U3fake[i] = U[i] + dt/2 * LU2[i];
    }
    double* LU3 = getLU(U3fake, x);  

    double* NU = malloc(U_SIZE * sizeof(double));
    for (size_t i = 0; i < U_SIZE; i++)
    {
        NU[i] = U[i] + dt * (1.0/6*LU[i] + 1.0/3*LU1[i] + 1.0/3*LU2[i] + 1.0/6*LU3[i]);
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

void print_matrix(double* src) {
    for (size_t i = 0; i < U_SIZE; i++)
    {
        printf("[%ld] = %f\n", i, src[i]);
    }
}

int main() {
    double* x = linspace(-1.0, 1.0, nNodos);

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

    //get first 3U 
    double* U1 = getNU(U, x);
    double* U2 = getNU(U1, x);
    double* U3 = getNU(U2, x);

    double* UNMenos3 = U;
    double* UNMenos2 = U1;
    double* UNMenos1 = U2;
    double* Un = U3;

    for (size_t i = 3400; i < 3600; i++)
    {
        printf("Un[%ld] = %f\n", i, Un[i]);
    }

    double* LUnMenos3 = getLU(UNMenos3, x);
    double* LUnMenos2 = getLU(UNMenos2, x);
    double* LUnMenos1 = getLU(UNMenos1, x);
    
    double tiempo = 3*dt;

    int m = 3;
    while(tMax > tiempo) {
        tiempo = tiempo + dt;
        double * LUn =  getLU(Un, x);

        double* UIM = malloc(U_SIZE * sizeof(double));
        for (size_t i = 0; i < U_SIZE; i++)
        {
            UIM[i] = Un[i] + dt*(55.0/24*LUn[i] - 59.0/24*LUnMenos1[i] + 37.0/24*LUnMenos2[i] -9.0/24*LUnMenos3[i]);
        }
        double* LUIM = getLU(UIM, x);
        
        double* UNmasUno = malloc(U_SIZE * sizeof(double));
        for (size_t i = 0; i < U_SIZE; i++)
        {
            UNmasUno[i] = Un[i] + dt*(9.0/24*LUIM[i] + 19.0/24*LUn[i] - 5.0/24*LUnMenos1[i] + 1.0/24*LUnMenos2[i]);
        }

        // double tol = 1e-15;
        double error = 1.0;
        int n = 1;

        while (error > 20*FLT_EPSILON)
        {
            double* LUNmasUno =  getLU(UNmasUno, x);
            double* UNmasUnoP = malloc(U_SIZE * sizeof(double));
            for (size_t i = 0; i < U_SIZE; i++)
            {
                UNmasUnoP[i] = Un[i] + dt*(9.0/24*LUNmasUno[i] + 19.0/24*LUn[i] - 5.0/24*LUnMenos1[i] + 1.0/24*LUnMenos2[i]);
            }

            
            double max = 0.0;
            for (size_t i = 0; i < U_SIZE; i++) {
                double resta = UNmasUno[i] - UNmasUnoP[i];
                double local_error = fabs(resta);
                if(local_error > max) {
                    max = local_error;
                }
            }

            error = max;

            for (size_t i = 0; i < U_SIZE; i++) {
                UNmasUno[i] = UNmasUnoP[i];
            }

            n = n + 1;

            free(LUNmasUno);
            free(UNmasUnoP);
        }

        for (size_t i = 0; i < U_SIZE; i++) {
            UNMenos3[i] = UNMenos2[i];
            UNMenos2[i] = UNMenos1[i];
            UNMenos1[i] = Un[i];
            Un[i] = UNmasUno[i];
        }

        LUnMenos3 = getLU(UNMenos3,x);
        LUnMenos2 = getLU(UNMenos2,x);
        LUnMenos1 = getLU(UNMenos1,x);
        
        // plot(x,Un(4,:),'-b.')
        // pause(0.01)
        sleep(0.01);
        m = m+1;

        free(UNmasUno);
        free(LUIM);
        free(UIM);
        free(LUn);
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
    free(U);
    free(U1);
    free(U2);
    free(U3);

    free(UNMenos3);
    free(UNMenos2);
    free(UNMenos1);
    free(Un);
    free(LUnMenos3);
    free(LUnMenos2);
    free(LUnMenos1);

    return 0;
}
