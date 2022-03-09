#include <iostream>
#include <stdio.h>
#include <windows.h>
#include <math.h>

using namespace std;


//  Zmienne globalne  
double TOL_BLEDU = 10e-10;
double TOL_FUN = 10e-10;  //dokladnosc, dokladnosc wartosci
int N_MAX = 50;

// funkcje  

double rownanie_funk_1(double x, double y, double z)
{
	return x * x + y * y + z * z - 2.0;
}

double rownanie_funk_2(double x, double y)
{
	return x * x + y * y - 1.0;
}

double rownanie_funk_3(double x, double y)
{
	return x * x - y;
}

void macierz_Jakobiego(double(&tab)[3][3], double x, double y, double z)
{
	tab[0][0] = 0.0;
	tab[0][1] = -2.0 * z;
	tab[0][2] = -4.0 * y * z;
	tab[1][0] = 0.0;
	tab[1][1] = -4.0 * x * z;
	tab[1][2] = 4.0 * x * z;
	tab[2][0] = 4.0 * x * z;
	tab[2][1] = 2.0 * x * (1.0 + 2.0 * y);
	tab[2][2] = 0.0;
}

double det(double x, double y, double z)
{
	return (-1.0 / ((4.0 * x * z) * (1.0 + 2.0 + y)));
}

void MacierzOdwrotana(double (&tab)[3][3], double x, double y, double z, double Mac_det)
{
	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			tab[i][j] *= Mac_det;
		}
	}
}

double norma(double tab[3])
{
	double max = fabs(tab[0]);
	for(int i = 1; i<3;i++)
	{
		if(fabs(tab[i])>max)
		{
			max = fabs(tab[i]);
		}
	}
	return max;
}

void wartBezwzgledna(double(&tab)[3], double(&wart)[3])
{
	for (int i = 0; i < 3; i++)
	{
		tab[i] = fabs(wart[i]);
	}
}




int main(void)
{
	cout << "METODY OBLICZENIOWE lab4 \n\n";

	double x = 0.4;
	double y = 0.4;
	double z = 0.9;


	int i = 0;
	double en = 0.1;
	double rezyduum = 0.1;
	
	double delta[3] = {0};
	double wektor[3] = {0};
	double Mac_det = 0.0;
	double M_Jakobiego[3][3] = {0};
	
	double ReBez[3] = {0};
	double EnBez[3] = {0};

	wektor[0] = rownanie_funk_1(x, y, z);
	wektor[1] = rownanie_funk_2(x, y);
	wektor[2] = rownanie_funk_3(x, y);

	cout << "|  NUMER|    X|     Y|     Z|        ESTYMATOR BLEDU|                REZYDUUM| \n";

	while(i++ < N_MAX && fabs(rezyduum) > TOL_FUN && fabs(en) > TOL_BLEDU)
	{
		Mac_det = det(x, y, z);
		macierz_Jakobiego(M_Jakobiego, x, y, z);
		MacierzOdwrotana(M_Jakobiego, x, y, z, Mac_det);
				
		delta[0] = 0;
		delta[1] = 0;
		delta[2] = 0;
	
		for(int j = 0; j < 3; j++)
		{
			for(int k = 0; k<3; k++)
			{
				delta[j] += M_Jakobiego[j][k] * wektor[k];
			}
		}
		
		wartBezwzgledna(ReBez, wektor);
		wartBezwzgledna(EnBez, delta);

		rezyduum = norma(ReBez);
		en = norma(EnBez);

		printf("|%d	|  %.5lf|  %.5lf|	%.5lf|    %.15lf|	    %.15lf|  \n", i, x, y, z, en, rezyduum);

		x -= delta[0];
		y -= delta[1];
		z -= delta[2];

		wektor[0] = rownanie_funk_1(x, y, z);
		wektor[1] = rownanie_funk_2(x, y);
		wektor[2] = rownanie_funk_3(x, y);
	}

	if (i > N_MAX)
    {
        cout << "\n # limit iteracji # \n";
    }
    else if (en < TOL_BLEDU)
    {
        cout << "\n # kryterium dokladnosci # \n";
    }
    else if (fabs(rezyduum) < TOL_FUN)
    {
        cout << "\n # kryterium wiarygodnosci  # \n";
    }

	return 0;
}


