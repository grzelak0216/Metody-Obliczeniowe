#include <iostream> 
using namespace std;

 //FLOAT 
float liczba_f = 1.0f;
float eps_f = liczba_f * 0.5f;
float tmp_f = liczba_f + 1.0f;
float epsf_maszyn;
int liczba_mantysy_f = 0;
 
 //DOUBLE 
double liczba_d = 1.0;
double eps_d = liczba_d * 0.5;
double tmp_d = liczba_d + 1.0;
double epsd_maszyn;
int liczba_mantysy_d = 0;
 
 void fun_FLOAT() 
 {
	while (tmp_f > liczba_f) 
	{
		epsf_maszyn= eps_f;
		eps_f = eps_f * 0.5f;
		tmp_f = liczba_f + eps_f;
		liczba_mantysy_f++;
	}
 }

void fun_DOUBLE()
{
	while (tmp_d > liczba_d) 
	{
		epsd_maszyn = eps_d;
		eps_d = eps_d * 0.5;
		tmp_d = liczba_d + eps_d;
		liczba_mantysy_d++; 
	}
}
 
int main() 
{
	fun_FLOAT(); 
	fun_DOUBLE(); 
	 
	cout << "\t\t FLOAT"<<endl;
	cout << "Liczba bitow santysy: " << liczba_mantysy_f << endl << "Epsilon maszynowy: " << epsf_maszyn << endl << endl;
	 
	cout << "\t\t DOUBLE"<<endl;
	cout << "Liczba bitow mantysy: " << liczba_mantysy_d << endl << "Epsilon maszynowy: " << epsd_maszyn << endl << endl;
	 
	system("pause"); 
	return 0;
}
