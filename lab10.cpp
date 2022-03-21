
#include<string>
#include<iostream>
#include<Windows.h>
#include<fstream>
#include<math.h>

using namespace std;

double rozwiazanie_analityczne(double t)
{
	return 1 - exp(-10.0 * (t + atan(t)));
}
double metoda_bezposrednia_Eulera(double h, double tmax, bool flag)
{
	double y = 0.0;
	double wd;
	double blad = 0.0;

	for (double t = 0.0; t < tmax; t += h)
	{
		wd = rozwiazanie_analityczne(t);
		y = y + h * (-((10.0 * t * t + 20.0) / (t * t + 1.0)) * (y - 1.0));
		wd = fabs(wd - y);

		if (wd > blad)
		{
			blad = wd;
		}
	}

	if (flag == true)
	{
		return y;
	}
	else
	{
		return blad;
	}
}

double metoda_posrednia_Eulera(double h, double tmax, bool flag)
{
	double y = 0.0, f;
	double wd;
	double blad = 0.0;

	for (double t = 0.0; t < tmax; t += h)
	{
		wd = rozwiazanie_analityczne(t);
		f = (10.0 * (t + h) * (t + h) + 20.0) / ((t + h) * (t + h) + 1.0);
		y = (y + h * f) / (1 + h * f);
		wd = fabs(wd - y);

		if (wd > blad)
		{
			blad = wd;
		}
	}	
	
	if (flag == true)
	{
		return y;
	}
	else
	{
		return blad;
	}
}

double metoda_trapezow(double h, double tmax, bool flag)
{
	double y = 0.0, f1, f2;
	double wd;
	double blad = 0.0;

	for (double t = 0.0;
		t < tmax;
		t += h) {
		wd = rozwiazanie_analityczne(t);
		y = ((-h / 2.0) * (((10.0 * t * t + 20.0) / (t * t + 1.0)) * (y - 1.0) - (10.0 * (t + h) * (t + h) + 20.0) / ((t + h) * (t + h) + 1.0)) + y) / (1.0 + (h / 2.0) * (10.0 * (t + h) * (t + h) + 20.0) / ((t + h) * (t + h) + 1.0));
		wd = fabs(wd - y);

		if (wd > blad)
		{
			blad = wd;
		}
	}
	
	if (flag == true)
	{
		return y;
	}
	else
	{
		return blad;
	}
}

int main()
{
	double krok = 0.1;
	const double t_max = 5.0;

	FILE* bledy = NULL;
	FILE* wyniki = NULL;
	FILE* BME_ns = NULL;

	bledy = fopen("mo10_bledy.txt", "w");
	if (bledy == NULL)
	{
		printf("Nie mozna otworzyc pliku");
		exit(1);
	}

	for (krok; krok > 1e-14; krok *= 0.5)
	{
		fprintf(bledy, "%lf\t%lf\t%lf\t%lf\n", log10(krok), log10(metoda_bezposrednia_Eulera(krok, t_max, false)), log10(metoda_posrednia_Eulera(krok, t_max, false)), log10(metoda_trapezow(krok, t_max, false)));
		printf("%.14lf \n", krok);
	}
	fclose(bledy);


	wyniki = fopen("mo10_wyniki.txt", "w");
	if (wyniki == NULL) {
		printf("Nie mozna otworzyc pliku");
		exit(1);
	}
	krok = 0.01;

	for (double t = 0; t < 5; t += 0.01)
	{
		fprintf(wyniki, "%lf\t%lf\t%lf\t%lf\t%lf\n", t, rozwiazanie_analityczne(t), metoda_bezposrednia_Eulera(krok, t, true), metoda_posrednia_Eulera(krok, t, true), metoda_trapezow(krok, t, true));
		printf("%.14lf \n", t);
	}
	fclose(wyniki);


	BME_ns = fopen("mo10_BME_niestabilne.txt", "w");
	if (BME_ns == NULL)
	{
		printf("Nie mozna otworzyc pliku");
		exit(1);
	}

	for (double t = 0; t < 5; t += 0.15)
	{
		fprintf(BME_ns, "%lf\t%lf\t%lf\n", t, rozwiazanie_analityczne(t), metoda_bezposrednia_Eulera(0.2, t, true));
		printf("%.14lf \n", t);
	}
	fclose(BME_ns);


	system("Pause");
	return 0;
}
