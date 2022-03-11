#include <iostream>
#include <cmath>
#include <fstream>
#include <cstdlib>
#include <iomanip>

using namespace std;

double **newMatrix(int n, int m)
{
    double **Matr;
    Matr = new double *[n];

    for (int i = 0; i < n; i++)
    {
        Matr[i] = new double[m];
    }

    return Matr;
}

void releaseMatrix(double **Matr, int n)
{
    for (int i = n; i > 0; i--)
    {
        delete[] Matr[i];
    }
    delete[] Matr;
}

void fillBasickMatrix(double **Matr)
{
    Matr[0][0] = 1.0;
    Matr[0][1] = -20.0;
    Matr[0][2] = 30.0;
    Matr[0][3] = -4.0;
    Matr[1][0] = 2.0;
    Matr[1][1] = -40.0;
    Matr[1][2] = -6.0;
    Matr[1][3] = 50.0;
    Matr[2][0] = 9.0;
    Matr[2][1] = -180.0;
    Matr[2][2] = 11.0;
    Matr[2][3] = -12.0;
    Matr[3][0] = -16.0;
    Matr[3][1] = 15.0;
    Matr[3][2] = -140.0;
    Matr[3][3] = 13.0;
}

void fillL(double **Matr)
{
    Matr[0][0] = 1.0;
    Matr[0][1] = 0.0;
    Matr[0][2] = 0.0;
    Matr[0][3] = 0.0;
    Matr[1][0] = 0.0;
    Matr[1][1] = 1.0;
    Matr[1][2] = 0.0;
    Matr[1][3] = 0.0;
    Matr[2][0] = 0.0;
    Matr[2][1] = 0.0;
    Matr[2][2] = 1.0;
    Matr[2][3] = 0.0;
    Matr[3][0] = 0.0;
    Matr[3][1] = 0.0;
    Matr[3][2] = 0.0;
    Matr[3][3] = 1.0;
}

void printMatrix(double **Matr, int n, int m)
{
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < m; j++)
        {
            cout << setw(5) << Matr[i][j] << " ";
        }
        cout << endl;
    }
    cout << endl;
}

void printMatrixL(double **Matr, int n, int m, int *r)
{
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < m; j++)
        {
            cout << setw(5) << Matr[r[i]][j] << " ";
        }
        cout << endl;
    }
    cout << endl;
}

double *newVector(int n)
{
    double *Vect;
    Vect = new double[n];

    return Vect;
}

void releaseVetor(double *Vect)
{
    delete[] Vect;
}

void fillBasickVector(double *Vect)
{
    Vect[0] = 35.0;
    Vect[1] = 104.0;
    Vect[2] = -366.0;
    Vect[3] = -354.0;
}

void printVector(double *Vect, int n)
{
    for (int i = 0; i < n; i++)
    {
        cout << setw(5) << Vect[i] << " ";
    }
    cout << endl;
}

void searchRow(double **Matr, int idi, int idm, int n)
{
    int j = n - 2;
    for (int i = 0; i < j; i++)
    {
        if (idi != i)
        {
            double k = Matr[idi][i];
            Matr[idi][i] = Matr[idm][i];
            Matr[idm][i] = k;
        }
    }
}

void GaussEl(double **Matr, double **L, int *r)
{
    int row;
    double tmp;

    for (int i = 0; i < 3; i++)
    {
        if (Matr[r[i]][r[i]] == 0.0)
        {
            int idmax;
            for (int j = i + 1; j < 4; j++)
            {
                double max = 0.0;
                if (fabs(Matr[r[j]][r[i]]) > max)
                {
                    max = Matr[r[j]][r[i]];
                    idmax = j;
                }
            }
            searchRow(L, i, idmax, 4);
            double k = r[idmax];
            r[idmax] = r[i];
            r[i] = k;
        }

        for (int j = i + 1; j < 4; j++)
        {
            tmp = Matr[r[j]][i] / Matr[r[i]][i];
            L[j][i] = tmp;
            for (int k = i; k < 4; k++)
            {
                Matr[r[j]][k] = Matr[r[j]][k] - Matr[r[i]][k] * tmp;
            }
        }
    }
}

void solution(double **Matr, double **L, double *vectorx, int *r, int n)
{
    cout << endl
         << "Macierz U:" << endl;
    printMatrixL(Matr, 4, 4, r);

    cout << endl
         << "Macierz U:" << endl;
    printMatrix(L, 4, 4);

    double *vectory = new double[n];

    for (int i = 0; i <= n; i++)
    {
        double sum = 0.0;

        for (int j = 0; j < i; j++)
        {
            sum = sum + L[i][j] * vectory[r[j]];
        }

        vectory[r[i]] = vectorx[r[i]] - sum;
    }

    for (int i = n; i >= 0; i--)
    {
        double sum = 0.0;

        for (int j = i + 1; j < 4; j++)
        {
            sum = sum + Matr[r[i]][j] * vectorx[r[j]];
        }

        vectorx[r[i]] = (1 / Matr[r[i]][i]) * (vectory[r[i]] - sum);
    }

    cout << "Rozwiazanie: " << endl
         << endl;

    for (int i = 0; i < 4; i++)
    {
        cout << "| " << setw(5) << vectorx[r[i]] << " |" << endl;
    }
}

int main()
{
    int rv[4] = {0, 1, 2, 3};

    double **Matrix = newMatrix(4, 4);
    fillBasickMatrix(Matrix);

    double *Vector = newVector(4);
    fillBasickVector(Vector);

    cout << "Macierz A: " << endl;
    printMatrix(Matrix, 4, 4);

    cout << "Wektor B: " << endl;
    printVector(Vector, 4);

    double **L = newMatrix(4, 4);
    fillL(L);

    GaussEl(Matrix, L, rv);

    solution(Matrix, L, Vector, rv, 3);

    releaseVetor(Vector);
    releaseMatrix(Matrix, 4);
    releaseMatrix(L, 4);

    system("Pause");
    return 0;
}