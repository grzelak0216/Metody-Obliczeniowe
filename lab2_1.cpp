#include <iostream> 
#include <string> 
#include <iomanip> 
#include <math.h> 
#include <fstream> 

using namespace std;

double fun_fx(double x)
{
	return (1 - exp(-x)) / x;
}

int main()
{
	ifstream indata;
	ofstream outdata;
	
	double x, f, f0, delx;
	string inlogx, inx, infx;
	
	indata.open("dane.txt"); 
	outdata.open("wynik.txt"); 
	if (!indata.is_open() || !outdata.is_open())
	{
		if (!outdata.is_open())
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
		f = fun_fx(x); 
		
		delx = log10(fabs(f - f0) / f0); 
		
		outdata << inlogx << " " << setprecision(21) << delx << scientific << endl; 
		
	}	
	
	indata.close(); 
	outdata.close();
	
	system("pause"); 
	return 0;
}
