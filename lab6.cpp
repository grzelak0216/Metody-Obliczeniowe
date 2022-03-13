#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <cmath>

using namespace std;

double *newVector(int n)
{
    double *Vect;
    Vect = new double[n];

    return Vect;
}

void vectorValue(double *v1, double *v2, double *v3, double *b)
{
    v1[0] = 1.0 / 2.0;
    v1[1] = 1.0 / 4.0;
    v1[2] = 1.0 / 6.0;
    v1[3] = 1.0 / 8.0;
    v1[4] = 1.0 / 10.0;
    v2[0] = 10;
    v2[1] = 20;
    v2[2] = 30;
    v2[3] = 30;
    v2[4] = 20;
    v2[5] = 10;
    v3[0] = 1.0 / 3.0;
    v3[1] = 1.0 / 5.0;
    v3[2] = 1.0 / 7.0;
    v3[3] = 1.0 / 9.0;
    v3[4] = 1.0 / 11.0;
    b[0] = 31;
    b[1] = 165.0 / 4.0;
    b[2] = 917.0 / 30.0;
    b[3] = 851.0 / 28.0;
    b[4] = 3637.0 / 90.0;
    b[5] = 332.0 / 11.0;
}

void printV(double *v, int n)
{
    for (int i = 0; i < n; i++)
    {
        printf("[ %5.5lf ]\n", v[i]);
    }
    cout << endl;
}

void removeVector(double *v1, double *v2, double *v3, double *v4, double *v5)
{
    delete[] v1;
    delete[] v2;
    delete[] v3;
    delete[] v4;
    delete[] v5;
}

void thomasMatrix(double *L, double *D, double *U)
{
    for (int i = 0; i < 5; i++)
    {
        D[i + 1] -= L[i] * (U[i] / D[i]);
        L[i] = L[i] / D[i];
    }
}

void thomasVector(double *L, double *b)
{
    for (int i = 1; i < 6; i++)
    {
        b[i] -= L[i - 1] * b[i - 1];
    }

    for (int i = 0; i < 5; i++)
    {
        L[i] = 0;
    }
}

void solution(double *L, double *D, double *U, double *b, double *X)
{
    thomasMatrix(L, D, U);
    thomasVector(L, b);

    X[5] = b[5] / D[5];
    for (int i = 4; i >= 0; i--)
    {
        X[i] = (b[i] - U[i] * X[i + 1]) / D[i];
    }
}

int main()
{
    double *vectorL = newVector(5);
    double *vectorD = newVector(6);
    double *vectorU = newVector(5);
    double *vectorb = newVector(6);
    double *vectorX = newVector(6);

    vectorValue(vectorU, vectorD, vectorL, vectorb);

    cout << "Wektor L: " << endl;
    printV(vectorL, 5);
    cout << "Wektor D: " << endl;
    printV(vectorD, 6);
    cout << "Wektor U: " << endl;
    printV(vectorU, 5);
    cout << "Wektor b: " << endl;
    printV(vectorb, 6);

    solution(vectorL, vectorD, vectorU, vectorb, vectorX);

    cout << endl
         << endl
         << "Modyfikacja wektorow" << endl;

    cout << "Wektor L: " << endl;
    printV(vectorL, 5);
    cout << "Wektor D: " << endl;
    printV(vectorD, 6);
    cout << "Wektor U: " << endl;
    printV(vectorU, 5);
    cout << "Wektor b: " << endl;
    printV(vectorb, 6);
    cout << "Wynik Wektor x: " << endl;
    printV(vectorX, 6);

    removeVector(vectorL, vectorD, vectorU, vectorb, vectorX);

    system("Pause");
    return 0;
}