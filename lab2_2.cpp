#include <iostream> 
#include <string> 
#include <iomanip> 
#include <math.h> 
#include <fstream> 

using namespace std;

double fun_fx(double x) 
{
	return (1- exp(-x)) / x;
}

double fun_Tayolra(double x) 
{
	double result = 1, before = 1;
	if (x < 0.2) 
	{
		for (int i = 1; i < 30; i++)
		{
			before *= (-1.0 * x)/(double(i)+1); 
			result += before;
		} 
		
		return result;		
	}
	else
	{
		return fun_fx(x);	
	} 
}

int main()
{
	ifstream indata;
	ofstream outdata2;
	
	double x, f, f0, delx;
	string inlogx, inx, infx;
	
	indata.open("dane.txt"); 
	outdata2.open("wynik2.txt"); 
	if (!indata.is_open() || !outdata2.is_open())
	{
		if (!outdata2.is_open())
		{
			cout << "Blad otwarcia pliku do zapisu" << endl;
		}
		else
		{
			cout << "Blad otwarcia pliku z danymi wejsciowymi"<<endl;
		} 
	}
	
	while (!indata.eof()) 
	{ 
		indata >> inlogx >> inx >> infx;
		x = std::stod(inx); 
		f0 = std::stod(infx);
		f = fun_Tayolra(x); 
		
		delx = log10(fabs(f - f0) / f0); 
		
		outdata2 << inlogx << " " << setprecision(21) << delx << scientific << endl; 
		
	}	
	
	indata.close(); 
	outdata2.close();
	
	system("pause"); 
	return 0;
}
