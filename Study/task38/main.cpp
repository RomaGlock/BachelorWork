#include <iostream>
#include <cmath>
#include <fstream>
using namespace std;

const double Pi = 3.141592653589;
const double T = 1e-1;
const int M = 30;
const int N = 30;
const int K = 65;

double U[M][N][2 * K];
double hx = 1.0 / (M - 1);
double hy = 2.0 / (N - 1);
double ht = T / (K - 1);

double F(int, int, int, char);
void boundary_conditions(int);
void printU(void);

int main(void)
{
    double A1 = ht / hx / hx / 2, C1 = A1, B1 = 1 + 2 * A1;
    double A2 = ht / hy / hy / 2, C2 = A2, B2 = 1 + 2 * A2;
    double dm[M], sm[M], dn[N], sn[N];
    dm[0] = 1; sm[0] = 0;
    dn[0] = 0; sn[0] = 0;

    for (int m = 0; m < M; m++)
        for (int n = 0; n < N; n++)
        U[m][n][0] = sin(2 * Pi * hx * m) * sin(Pi * 0.25 * hy * n);

    for (int k = 0; k < 2 * K - 3; k += 2)
    {
        for (int n = 1; n < N - 1; n++)
        {
            for (int m = 0; m < M - 1; m++)
            {
                dm[m + 1] = C1 / (B1 - A1*dm[m]);
                sm[m + 1] = (F(m, n, k, '1') + A1*sm[m]) / (B1 - A1*dm[m]);
            }

            U[M - 1][n][k + 1] = sm[M - 1] / (1 - dm[M - 1]);

            for (int m = M - 1; m > 1; m--)
                U[m - 1][n][k + 1] = dm[m] * U[m][n][k + 1] + sm[m];
        }

        boundary_conditions(k);

        for (int m = 1; m < M - 1; m++)
        {
            for (int n = 0; n < N - 1; n++)
            {
                dn[n + 1] = C2 / (B2 - A2*dn[n]);
                sn[n + 1] = (F(m, n, k, '2') + A2*sn[n]) / (B2 - A2*dn[n]);
            }

            U[m][N - 1][k + 2] = 0;

            for (int n = N - 1; n > 1; n--)
                U[m][n - 1][k + 2] = dn[n] * U[m][n][k + 2] + sn[n];
        }

        boundary_conditions(k+1);
    }

    printU();

}

void printU(void)
{
    ofstream foutMeshConfig("mesh.config");
    ofstream foutU0("U0.dat");
    ofstream foutU32("U32.dat");
    ofstream foutU64("U64.dat");
    ofstream foutU128("U128.dat");



    foutMeshConfig << "Hx = " << fixed << hx << " N = " << M << endl;
    foutMeshConfig << "Hy = " << fixed << hy << " M = " << N << endl;

    for (int m = 0; m < M; m++)
    {
        for (int n = 0; n < N; n++)
        {
            foutU0      << U[m][n][0]   << " ";
            foutU32     << U[m][n][32]  << " ";
            foutU64     << U[m][n][64]  << " ";
            foutU128    << U[m][n][128] << " ";
        }
        foutU0      << endl;
        foutU32     << endl;
        foutU64     << endl;
        foutU128    << endl;
    }
}

double F(int m, int n, int k, char s)
{
    double a = 0.0;

    switch (s)
    {

    case '1':
        a = (U[m][n - 1][k] + U[m][n + 1][k])*ht / hy / hy / 2 + (1 - ht / hy /
        hy)*U[m][n][k];
        break;

    case '2':
        a = (U[m - 1][n][k + 1] + U[m + 1][n][k + 1])*ht / hx / hx / 2 + (1 - ht /
        hx / hx)*U[m][n][k + 1];
        break;

    }

    return a;
}

void boundary_conditions(int k)
{
    for (int n = 0; n < N; n++)
    {
        U[0][n][k + 1] = 0;
        U[M - 1][n][k + 1] = 0;
    }

    for (int m = 0; m < M; m++)
    {
        U[m][0][k + 1] = 0;
        U[m][N - 1][k + 1] = U[m][N - 2][k + 1];
    }
}
