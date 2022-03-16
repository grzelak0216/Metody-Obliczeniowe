#include <iostream>
#include <math.h>
#include <iomanip>
#include <stdlib.h>

using namespace std;

double EPS = 1e-8;
const int loop = 50;

double** newMatrix(int n, int m)
{
    double **Matr;
    Matr = new double *[n];

    for (int i = 0; i < n; i++)
    {
        Matr[i] = new double[m];
    }

    return Matr;
}

double *newVector(int n)
{
    double *Vect;
    Vect = new double[n];
    
    for (int i = 0; i < n; i++)
    {
    	Vect[i] = 0.0;
	}

    return Vect;
}

void fillBasickMatrix(double **Matr)
{
    Matr[0][0] = 100.0;
    Matr[0][1] = -1.0;
    Matr[0][2] = 2.0;
    Matr[0][3] = -3.0;
    Matr[1][0] = 1.0;
    Matr[1][1] = 200.0;
    Matr[1][2] = -4.0;
    Matr[1][3] = 5.0;
    Matr[2][0] = -2.0;
    Matr[2][1] = 4.0;
    Matr[2][2] = 300.0;
    Matr[2][3] = -6.0;
    Matr[3][0] = 3.0;
    Matr[3][1] = -5.0;
    Matr[3][2] = 6.0;
    Matr[3][3] = 400.0;
}

void fillVectorB(double *Vect)
{
    Vect[0] = 116.0;
    Vect[1] = -226.0;
    Vect[2] = 912.0;
    Vect[3] = -1174.0;
}

void fillVectorX(double *Vect)
{
    Vect[0] = 2;
    Vect[1] = 2;
    Vect[2] = 2;
    Vect[3] = 2;
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

void printVector(double *Vect, int n)
{
    for (int i = 0; i < n; i++)
    {
        cout << setw(5) << Vect[i] << " ";
    }
    cout << endl;
}

void releaseMatrix(double **Matr, int n)
{
    for (int i = n; i > 0; i--)
    {
        delete[] Matr[i];
    }
    delete[] Matr;
}

void releaseVetor(double *Vect)
{
    delete[] Vect;
}

double estymator(double* x, double* x_nowe)
{
	x[0] = fabs(x[0] - x_nowe[0]);
	x[1] = fabs(x[1] - x_nowe[1]);
	x[2] = fabs(x[2] - x_nowe[2]);
	x[3] = fabs(x[3] - x_nowe[3]);
	
	double max = x[0];
	
	for (int i = 1; i < 4; i++)
	{
		if(x[i] > max)
		{
			max = x[i];
		}
	}
	return max;
}

double residuum(double** m, double* b, double* x_nowe)
{
	double Ax[4];
	
	Ax[0] = fabs((m[0][0] * x_nowe[0] + m[0][1] * x_nowe[1] + m[0][2] * x_nowe[2] + m[0][3] * x_nowe[3]) - b[0]);
	Ax[1] = fabs((m[1][0] * x_nowe[0] + m[1][1] * x_nowe[1] + m[1][2] * x_nowe[2] + m[1][3] * x_nowe[3]) - b[1]);
	Ax[2] = fabs((m[2][0] * x_nowe[0] + m[2][1] * x_nowe[1] + m[2][2] * x_nowe[2] + m[2][3] * x_nowe[3]) - b[2]);
	Ax[3] = fabs((m[3][0] * x_nowe[0] + m[3][1] * x_nowe[1] + m[3][2] * x_nowe[2] + m[3][3] * x_nowe[3]) - b[3]);
	
	double max = Ax[0];
	
	for (int i = 1; i < 4; i++)
	{
		if (Ax[i] > max)
		{
			max = Ax[i];
		}
	}
	
	return max;
}


void jacobi(double** m, double* b, double* x, const int n)
{
	double s, est = 1.0, res = 1.0;
	int ii = 0;
	double* x2 = new double[n];
	
	cout << "| lp.| x[1] |  x[2] |  x[3] |  x[4] | estymator | resyduum|\n";
	
	while(1)
	{
		for (int i = 0; i < n; i++)
		{
			s = 0.0;
			for (int j = 0; j < n; j++)
			{
				if (i != j)
				{
					s += m[i][j] * x[j];
				}
			}
			x2[i] = (b[i] - s) / m[i][i];
		}
		
		est = estymator(x, x2);
		res = residuum(m, b, x2);
		
		for (int i = 0; i < n; i++)
		{
			x[i] = x2[i];
		}
		
		ii++;
		
		printf("|%3.0d| %2.10e| %2.10e| %2.10e| %2.10e| %2.10e| %2.10e|\n", ii, x[0], x[1], x[2], x[3], est, res);
		
		if ((est < EPS && res < EPS) || ii > loop)
		{
			if (ii > loop)
			{
				cout << endl << "Przekroczono ilosc iteracji" << endl;
			}
			else
			{
				cout << endl << "Przekroczono kryterium dokladnosci" << endl;
			}
			break;
		}
	}
	delete[] x2;
}


void gauss_seidel(double** m, double* b, double* x, const int n)
{
	double s, est = 1.0, res = 1.0;
	int ii = 0;
	double* x2 = new double[n];
	
	cout << "| lp.| x[1] |  x[2] |  x[3] |  x[4] | estymator | resyduum|\n";
	
	while(1)
	{
		for (int i = 0; i < n; i++)
		{
			s = 0.0;
			for (int j = 0; j < n; j++)
			{
				if (i != j)
				{
					s += m[i][j] * x[j];
				}
			}
			x2[i] = x[i];
			x[i] = (b[i] - s) / m[i][i];
		}
		
		est = estymator(x2, x);
		res = residuum(m, b, x);
		
		ii++;
		
		printf("|%3.0d| %2.10e| %2.10e| %2.10e| %2.10e| %2.10e| %2.10e|\n", ii, x[0], x[1], x[2], x[3], est, res);
		
		if ((est < EPS && res < EPS) || ii > loop)
		{
			if (ii > loop)
			{
				cout << endl << "Przekroczono ilosc iteracji" << endl;
			}
			else
			{
				cout << endl << "Przekroczono kryterium dokladnosci" << endl;
			}
			break;
		}
	}
	delete[] x2;
}

void sor(double** m, double* b, double* x, double omega, const int n)
{
	double s, est = 1.0, res = 1.0;
	int ii = 0;
	double* x2 = new double[n];
	double* x3 = new double[n];
	
	cout << "| lp.| x[1] |  x[2] |  x[3] |  x[4] | estymator | resyduum|\n";
	
	while(1)
	{
		for (int i = 0; i < n; i++)
		{
			s = 0.0;
			for (int j = 0; j < n; j++)
			{
				if (i != j)
				{
					s += m[i][j] * x[j];
				}
			}
			x3[i] = x[i];
			x2[i] = (1.0 - omega) * x[i] + (omega / m[i][i]) * (b[i] - s);
			x[i] = x2[i];
		}
		
		est = estymator(x, x3);
		res = residuum(m, b, x2);
		
		for (int i = 0; i < n; i++)
		{
			x[i] = x2[i];
		}
		
		ii++;
		
		printf("|%3.0d| %2.10e| %2.10e| %2.10e| %2.10e| %2.10e| %2.10e|\n", ii, x[0], x[1], x[2], x[3], est, res);
		
		if ((est < EPS && res < EPS) || ii > loop)
		{
			if (ii > loop)
			{
				cout << endl << "Przekroczono ilosc iteracji" << endl;
			}
			else
			{
				cout << endl << "Przekroczono kryterium dokladnosci" << endl;
			}
			break;
		}
	}
	delete[] x2;
}
	
int menu_start()
{
	int ms;
		
	cout << "Rozwiazanie ukladow liniowych metodami iteracyjnymi" << endl << endl
		<<"Wybierz interesujaca Cie metode: " << endl
		<<"1. Metoda Jacobiego" << endl
		<<"2. Metoda Gaussa-Seidela" << endl
		<<"3. Metoda SOR" << endl 
		<<"4. Wyjdz z programu" << endl << endl;
		
	cin >> ms;
		
	while (ms < 1 || ms > 4)
	{
		cout << "zla wartosc" << endl;
		cin >> ms;
	}
		
	return ms;								
}

int main()
{
	double** Matrix = newMatrix(4, 4);
	double* VectorB = newVector(4);
	double* VectorX = newVector(4);
	
	fillBasickMatrix(Matrix);
	fillVectorB(VectorB);
	fillVectorX(VectorX);
	
	bool exit = false;
	
	while (!exit)
	{
		cout << "Macierz A: " << endl;
		printMatrix(Matrix, 4, 4);
		
		cout << "Wektor B: " << endl;
		printVector(VectorB, 4);
		
		cout << "Wektor X: " << endl;
		printVector(VectorX, 4);
		
		switch (menu_start())
		{
		case 1:
			cout << endl << "Wybrano: Metoda Jacobiego " << endl << endl;
			jacobi(Matrix, VectorB, VectorX, 4);
			system("pause");
			break;
		
		case 2:
			cout << endl << "Wybrano: Metoda Gaussa-Seidela " << endl << endl;
			gauss_seidel(Matrix, VectorB, VectorX, 4);
			system("pause");
			break;
			
		case 3:
			cout << endl << "Wybrano: Metoda SOR " << endl << endl;
			sor(Matrix, VectorB, VectorX, 0.5, 4);
			system("pause");
			break;
			
		case 4:
			cout << endl << "Koniec programu " << endl << endl;
			exit = true;
			system("pause");
			break;
		}
		
		if (!exit)
		{
			fillVectorX(VectorX);
			system("cls");
		}
	}
	releaseVetor(VectorB);
	releaseVetor(VectorX);
	releaseMatrix(Matrix, 4);
	
	system("pause");
	return 0;
}	
