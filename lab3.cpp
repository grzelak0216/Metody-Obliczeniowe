#include <iostream>
#include <stdio.h>
#include <windows.h>
#include <math.h>

using namespace std;

double TOL_BLEDU = 10e-10, TOL_FUN = 10e-10;
int N_MAX = 50;

double rownanie_funk_A(double x)
{
	double wynik = sin(x / 4.0);
	wynik *= wynik;
	return wynik - x;
}

double fi_funk_A(double x)
{
	double wynik = sin(x / 4.0);
	return wynik * wynik;
}

double pochodna_funk_A(double x)
{
	double wynik = 0.5 * sin(x / 4.0) * cos(x / 4.0) - 1;
	return wynik;
}


double rownanie_funk_B(double x)
{
	double wynik = tan(2.0 * x) - x - 1.0;
	return wynik;
}

double fi_funk_B(double x)
{
	double wynik = tan(2.0 * x) - 1.0;
	return wynik;
}

double pochodna_funk_B(double x)
{
	double wynik = cos(2.0 * x);
	wynik *= wynik;
	return (0.5 / wynik);
}




void picardA()
{
	cout << "\n---------------------------Metoda Picarda---------------------------\n\n";
	
	int i = 0;
	double x = 1.0;
	double x0 = 0;
	
	cout << "|numer  |przyblizenie       |estymator bledu     | rezyduum      |\n";
	
	while(i++ <= N_MAX && fabs(x - x0) > TOL_BLEDU && fabs(rownanie_funk_A(x)) > TOL_FUN)
	{
		x0 = x;
		x = fi_funk_A(x);
		
		printf("|%d  | %.15lf|  %.15lf| %.15lf| \n", i, x, fabs(x-x0), fabs(rownanie_funk_A(x)));
	}
	
	cout << "\n-------------------------------koniec-------------------------------\n\n";
	
	
	if (i++ > N_MAX)
	{
		cout << "\n # limit iteracji # \n";
	} 
	else if (fabs(x - x0) < TOL_BLEDU)
	{
		cout << "\n # kryterium dokladnosci # \n";
	}
	else if (fabs(rownanie_funk_A(x)) > TOL_FUN)
	{
		cout << "\n # kryterium wiarygodnosci # \n";
	}
	
	getchar();
}


void bisekcjaA()
{
	cout << "\n---------------------------Metoda Bisekcja---------------------------\n\n";
	
	int i = 0;
	double a = -0.2, b = 1.0, c;
	
	cout << "|numer  |przyblizenie       |estymator bledu     | rezyduum      |\n";
	
	while(i++ <= N_MAX && rownanie_funk_A(a) != 0 && rownanie_funk_A(b) != 0 && fabs((b - a) / 2.) >= TOL_BLEDU && fabs(rownanie_funk_A(c)) > TOL_FUN)
	{
		c = (a + b) / 2;
		
		if((rownanie_funk_A(a) < 0 && rownanie_funk_A(c) > 0) || (rownanie_funk_A(a) > 0 && rownanie_funk_A(c) < 0))
		{
			b = c;
		}
		else
		{
			if((rownanie_funk_A(c) < 0 && rownanie_funk_A(b) > 0) || (rownanie_funk_A(c) > 0 && rownanie_funk_A(b) < 0))
			{
				a = c;
			}
			else
			{
				printf("Brak wyniku");
				break;
			}
		}
		
		printf("|%d  | %.15lf|  %.15lf| %.15lf| \n", i, c, fabs((b - a) / 2.), fabs(rownanie_funk_A(c)));
	}
	
	cout << "\n-------------------------------koniec-------------------------------\n\n";
	
	
	if (i > N_MAX)
	{
		cout << "\n # limit iteracji # \n";
	} 
	else if (fabs((b - a) / 2.) < TOL_BLEDU)
	{
		cout << "\n # kryterium dokladnosci # \n";
	}
	else if (fabs(rownanie_funk_A(c)) > TOL_FUN)
	{
		cout << "\n # kryterium wiarygodnosci # \n";
	}
	
	getchar();
}


void newtonA()
{
	cout << "\n---------------------------Metoda Newtona---------------------------\n\n";
	
	int i = 0;
	double x = 0.1;
	double x0 = 0.0;
	
	cout << "|numer  |przyblizenie       |estymator bledu     | rezyduum      |\n";
	
	while(i <= N_MAX && fabs(x - x0) > TOL_BLEDU && fabs(rownanie_funk_A(x)) > TOL_FUN)
	{
		x0 = x;
		x = x0 - (rownanie_funk_A(x0) / pochodna_funk_A(x0));
		
		printf("|%d  | %.15lf|  %.15lf| %.15lf| \n", i, x, fabs(x-x0), fabs(rownanie_funk_A(x)));
		i++;
	}
	
	cout << "\n-------------------------------koniec-------------------------------\n\n";
	
	
	if (i++ > N_MAX)
	{
		cout << "\n # limit iteracji # \n";
	} 
	else if (fabs(x - x0) < TOL_BLEDU)
	{
		cout << "\n # kryterium dokladnosci # \n";
	}
	else if (fabs(rownanie_funk_A(x)) > TOL_FUN)
	{
		cout << "\n # kryterium wiarygodnosci # \n";
	}
	
	getchar();
}


void sieczneA()
{
	cout << "\n---------------------------Metoda sieczne---------------------------\n\n";
	
	int i = 0;
	double x1;
	double x3 = 0.2, x2 = 1.0;
	
	cout << "|numer  |przyblizenie       |estymator bledu     | rezyduum      |\n";
	
	while(i++ <= N_MAX && fabs(x3 - x2) > TOL_BLEDU && fabs(rownanie_funk_A(x3)) > TOL_FUN)
	{
		x1 = x2;
		x2 = x3;
		x3 = x2 - (rownanie_funk_A(x2) * (x2 - x1)) / (rownanie_funk_A(x2) - rownanie_funk_A(x1));
		
		printf("|%d  | %.15lf|  %.15lf| %.15lf| \n", i, x3, fabs(x3-x2), fabs(rownanie_funk_A(x3)));
	}
	
	cout << "\n-------------------------------koniec-------------------------------\n\n";
	
	
	if (i++ > N_MAX)
	{
		cout << "\n # limit iteracji # \n";
	} 
	else if (fabs(x3 - x2) < TOL_BLEDU)
	{
		cout << "\n # kryterium dokladnosci # \n";
	}
	else if (fabs(rownanie_funk_A(x3)) > TOL_FUN)
	{
		cout << "\n # kryterium wiarygodnosci # \n";
	}
	
	getchar();
}


void picardB()
{
	cout << "\n---------------------------Metoda Picarda---------------------------\n\n";
	
	int i = 0;
	double x = 1.0;
	double x0 = 0;
	
	cout << "|numer  |przyblizenie       |estymator bledu     | rezyduum      |\n";
	
	while(i++ <= N_MAX && fabs(x - x0) > TOL_BLEDU && fabs(rownanie_funk_B(x)) > TOL_FUN)
	{
		x0 = x;
		x = fi_funk_B(x);
		
		printf("|%d  | %.15lf|  %.15lf| %.15lf| \n", i, x, fabs(x-x0), fabs(rownanie_funk_B(x)));
	}
	
	cout << "\n-------------------------------koniec-------------------------------\n\n";
	
	
	if (i++ > N_MAX)
	{
		cout << "\n # limit iteracji # \n";
	} 
	else if (fabs(x - x0) < TOL_BLEDU)
	{
		cout << "\n # kryterium dokladnosci # \n";
	}
	else if (fabs(rownanie_funk_B(x)) > TOL_FUN)
	{
		cout << "\n # kryterium wiarygodnosci # \n";
	}
	
	getchar();
}


void bisekcjaB()
{
	cout << "\n---------------------------Metoda Bisekcja---------------------------\n\n";
	
	int i = 0;
	double a = -0.2, b = 1.0, c;
	
	cout << "|numer  |przyblizenie       |estymator bledu     | rezyduum      |\n";
	
	while(i++ <= N_MAX && rownanie_funk_B(a) != 0 && rownanie_funk_B(b) != 0 && fabs((b - a) / 2.) >= TOL_BLEDU && fabs(rownanie_funk_B(c)) > TOL_FUN)
	{
		c = (a + b) / 2.;
		
		if((rownanie_funk_B(a) < 0 && rownanie_funk_B(c) > 0) || (rownanie_funk_B(a) > 0 && rownanie_funk_B(c) < 0))
		{
			
		}
		else
		{
			if((rownanie_funk_B(c) < 0 && rownanie_funk_B(b) > 0) || (rownanie_funk_B(c) > 0 && rownanie_funk_B(b) < 0))
			{
				a = c;
			}
			else
			{
				printf("Brak wyniku");
				break;
			}
		}
		
		printf("|%d  | %.15lf|  %.15lf| %.15lf| \n", i, c, fabs((b - a) / 2.), fabs(rownanie_funk_B(c)));
	}
	
	cout << "\n-------------------------------koniec-------------------------------\n\n";
	
	
	if (i > N_MAX)
	{
		cout << "\n # limit iteracji # \n";
	} 
	else if (fabs((b - a) / 2.) < TOL_BLEDU)
	{
		cout << "\n # kryterium dokladnosci # \n";
	}
	else if (fabs(rownanie_funk_B(c)) > TOL_FUN)
	{
		cout << "\n # kryterium wiarygodnosci # \n";
	}
	
	getchar();
}


void newtonB()
{
	cout << "\n---------------------------Metoda Newtona---------------------------\n\n";
	
	int i = 0;
	double x = 0.1;
	double x0 = 0.0;
	
	cout << "|numer  |przyblizenie       |estymator bledu     | rezyduum      |\n";
	
	while(i <= N_MAX && fabs(x - x0) > TOL_BLEDU && fabs(rownanie_funk_B(x)) > TOL_FUN)
	{
		x0 = x;
		x = x0 - (rownanie_funk_B(x0) / pochodna_funk_B(x0));
		
		printf("|%d  | %.15lf|  %.15lf| %.15lf| \n", i, x, fabs(x-x0), fabs(rownanie_funk_B(x)));
		i++;
	}
	
	cout << "\n-------------------------------koniec-------------------------------\n\n";
	
	
	if (i++ > N_MAX)
	{
		cout << "\n # limit iteracji # \n";
	} 
	else if (fabs(x - x0) < TOL_BLEDU)
	{
		cout << "\n # kryterium dokladnosci # \n";
	}
	else if (fabs(rownanie_funk_B(x)) > TOL_FUN)
	{
		cout << "\n # kryterium wiarygodnosci # \n";
	}
	
	getchar();
}


void sieczneB()
{
	cout << "\n---------------------------Metoda sieczne---------------------------\n\n";
	
	int i = 0;
	double x1;
	double x3 = 0.2, x2 = 1.0;
	
	cout << "|numer  |przyblizenie       |estymator bledu     | rezyduum      |\n";
	
	while(i++ <= N_MAX && fabs(x3 - x2) > TOL_BLEDU && fabs(rownanie_funk_B(x3)) > TOL_FUN)
	{
		x1 = x2;
		x2 = x3;
		x3 = x2 - (rownanie_funk_B(x2) * (x2 - x1)) / (rownanie_funk_B(x2) - rownanie_funk_B(x1));
		
		printf("|%d  | %.15lf|  %.15lf| %.15lf| \n", i, x3, fabs(x3-x2), fabs(rownanie_funk_B(x3)));
	}
	
	cout << "\n-------------------------------koniec-------------------------------\n\n";
	
	
	if (i++ > N_MAX)
	{
		cout << "\n # limit iteracji # \n";
	} 
	else if (fabs(x3 - x2) < TOL_BLEDU)
	{
		cout << "\n # kryterium dokladnosci # \n";
	}
	else if (fabs(rownanie_funk_B(x3)) > TOL_FUN)
	{
		cout << "\n # kryterium wiarygodnosci # \n";
	}
	
	getchar();
}


int main(void)
{
	cout << "1. sin2(x/4)-x \n\n";
	
	picardA();
	bisekcjaA();
	newtonA();
	sieczneA();
	
	getchar();
	
	cout << "\n\n2. tan(2x)-x-1 \n\n";
	
	picardB();
	bisekcjaB();
	newtonB();
	sieczneB();
}
