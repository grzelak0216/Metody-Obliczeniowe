#include <iostream>
#include <cmath>
#include <fstream>
#include "calerf.h"
#include "calerf.cpp"
using namespace std;

const double TOL = 1.0e-16;
const int N_MAX = 1000;

typedef void (*Method)(double **, double *, double *, int);

double estimator(double *old_x, double *new_x, int n)
{
	double *x = new double[n];

	for (int i = 0; i < n; i++)
	{
		x[i] = fabs(old_x[i] - new_x[i]);
	}

	double maxRes = x[0];

	for (int i = 1; i < n; i++)
	{
		if (x[i] > maxRes)
		{
			maxRes = x[i];
		}
	}

	delete[] x;
	return maxRes;
}

double residuum(double **A, double *b, double *new_x, int n)
{
	double *x = new double[n];

	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			x[i] += A[i][j] * new_x[i];
		}
		x[i] = fabs(x[i] - b[i]);
	}

	double maxRes = x[0];

	for (int i = 1; i < n; i++)
	{
		if (x[i] > maxRes)
		{
			maxRes = x[i];
		}
	}

	delete[] x;
	return maxRes;
}

void SOR(double **A, double *b, double *x, int n)
{
	double *x_r = new double[n];
	double *l_x = new double[n];
	double summation;
	double est;
	double res;
	double omega = 0.5;

	int iter = 0;
	do
	{
		for (int i = 0; i < n; i++)
		{
			summation = 0.0;
			for (int j = 0; j < n; j++)
			{
				if (j != i)
				{
					summation += A[i][j] * x[j];
				}
			}

			l_x[i] = x[i];
			x_r[i] = (1.0 - omega) * x[i] + (omega / A[i][i]) * (b[i] - summation);
			x[i] = x_r[i];
		}

		est = estimator(l_x, x_r, n);
		res = residuum(A, b, x_r, n);
	} while (iter < N_MAX && est > TOL && res > TOL);
}

void thomas(double **A, double *b, double *x, int n)
{
	double *d = new double[n];
	double *u = new double[n - 1];
	double *l = new double[n - 1];

	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			if (i == j)
			{
				d[i] = A[i][j];
				if (i == 0)
				{
					u[i] = A[i][j + 1];
				}
				else if (i == n - 1)
				{
					l[i - 1] = A[i][j - 1];
				}
				else
				{
					u[i] = A[i][j + 1];
					l[i - 1] = A[i][j - 1];
				}
			}
		}
	}

	double *tempb = new double[n];
	double *tempd = new double[n];

	tempd[0] = d[0];
	tempb[0] = b[0];

	for (int i = 1; i < n; i++)
	{
		tempd[i] = d[i] - l[i - 1] * (u[i - 1] / tempd[i - 1]);
	}

	for (int i = 1; i < n; i++)
	{
		tempb[i] = b[i] - l[i - 1] * tempb[i - 1] / tempd[i - 1];
	}

	x[n - 1] = tempb[n - 1] / tempd[n - 1];

	for (int i = n - 2; i >= 0; i--)
	{
		x[i] = (tempb[i] - u[i] * x[i + 1]) / tempd[i];
	}

	delete[] d;
	delete[] u;
	delete[] l;
	delete[] tempb;
	delete[] tempd;
}

void KMB(double *U, double lambda, int n)
{
	double *new_U = new double[n];
	for (int i = 1; i < n - 1; i++)
	{
		new_U[i] = lambda * U[i - 1] + (1.0 - 2.0 * lambda) * U[i] + lambda * U[i + 1];
	}
	new_U[0] = 1.0;
	new_U[n - 1] = 0.0;

	for (int i = 0; i < n; i++)
	{
		U[i] = new_U[i];
	}

	delete[] new_U;
}

void CrankNicolson_B(double *b, const double THETA, const double GAMMA, double *U, int n, double lambda)
{
	b[0] = -GAMMA;
	for (int i = 1; i < n - 1; i++)
	{
		b[i] = -(lambda * 0.5 * U[i - 1] + (1.0 - lambda) * U[i] + lambda * 0.5 * U[i + 1]);
	}
	b[n - 1] = -THETA;
}

void CrankNicolson_A(double **A, double lambda, double ALPHA, double BETA, double FI, double PSI, int n, double h)
{
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			if (i == j)
			{
				if (i == 0)
				{
					A[i][j] = -ALPHA / h + BETA;
					A[i][j + 1] = ALPHA / h;
				}
				else if (i == n - 1)
				{
					A[i][j] = -FI / h + PSI;
					A[i][j - 1] = FI / h;
				}
				else
				{
					A[i][j - 1] = lambda / 2.0;
					A[i][j] = -(1.0 + lambda);
					A[i][j + 1] = lambda / 2.0;
				}
			}
			else if (j != i + 1 && j != i - 1)
			{
				A[i][j] = 0.0;
			}
		}
	}
}

void CrankNicolson(Method method, double **A, double *b, double *U, int n)
{
	double *new_U = new double[n];
	for (int i = 0; i < n; i++)
	{
		new_U[i] = U[i];
	}

	method(A, b, new_U, n);

	for (int i = 0; i < n; i++)
	{
		U[i] = new_U[i];
	}

	delete[] new_U;
}

double *initU(double h, double a, double b, int n)
{
	double *U = new double[n];

	for (int i = 0; i < n; i++)
	{
		if (a + h * i < 0)
		{
			U[i] = 1.0;
		}
		else
		{
			U[i] = 0.0;
		}
	}
	return U;
}

void analyU(double *U, double h, double a, double b, int n, double t)
{
	for (int i = 0; i < n; i++)
	{
		U[i] = 0.5 * calerf::ERFCL((a + i * h) / (2.0 * sqrt(t)));
	}
}

double errorMax(double *U, double *analyU, int n)
{
	double maxErr = fabs(analyU[0] - U[0]);
	double temp;
	for (int i = 1; i < n; i++)
	{
		temp = fabs(analyU[i] - U[i]);
		if (maxErr < temp)
		{
			maxErr = temp;
		}
	}
	return maxErr;
}

int main()
{
	double tmax = 2.0;
	double D = 1.0;
	double a = -(6.0 * sqrt(D * tmax));
	double b = 6.0 * sqrt(D * tmax);
	double lambda_d = 0.4;	 //lambda for direct method
	double lambda_ind = 1.0; //lambda for indirect method

	double h;
	double dt_ind, dt_d; //dt (dt_ind - indirect, dt_d - direct)

	double ALPHA = 0.0;
	double BETA = 1.0;
	double FI = 0.0;
	double PSI = 1.0;
	double THETA = 0.0;
	double GAMMA = -1.0;

	int n;			//num of columns
	int k_ind, k_d; //num of rows (k_ind - indirect, k_d - direct)

	double H[10] = {0.5, 0.45, 0.40, 0.35, 0.30, 0.25, 0.20, 0.15, 0.10, 0.05}; 

    fstream file;
    file.open("errorTmax.txt", ios::out);
    file << "log10(h) \t log10(err) CN_Gaussa_Seidla \t log10(err) CN_THOMAS \t log10(err) KMB" << endl;
    for (int i = 0; i < 10; i++)
	{
		n = fabs(a - b) / H[i] + 1;
		//for KMB
		dt_d = (lambda_d * H[i] * H[i]) / D;
		k_d = (tmax / dt_d) + 2;
		//for Crank-Nicolson
		dt_ind = (lambda_ind * H[i] * H[i]) / D;
		k_ind = (tmax / dt_ind) + 2;
		//U init
		double *U = new double[n];
		//analyU
		double *aU = new double[n];
		analyU(aU, H[i], a, b, n, tmax);
		//init b
		double *B = new double[n];

		//==================== Matrix Crank-Nicolson ====================
		double **A = new double *[n];
		for (int i = 0; i < n; i++)
		{
			A[i] = new double[n];
		}
		CrankNicolson_A(A, lambda_ind, ALPHA, BETA, FI, PSI, n, H[i]);

		//==================== Crank-Nicolson with SOR ====================
		U = initU(H[i], a, b, n);
		for (int k_ = 1; k_ < k_ind; k_++)
		{
			CrankNicolson_B(B, THETA, GAMMA, U, n, lambda_ind);
			CrankNicolson(SOR, A, B, U, n);
		}
		file << scientific << log10(H[i]) << "  " << log10(errorMax(U, aU, n));

		//==================== Crank-Nicolson with Thomas ====================
		U = initU(H[i], a, b, n);
		for (int k_ = 1; k_ < k_ind; k_++)
		{
			CrankNicolson_B(B, THETA, GAMMA, U, n, lambda_ind);
			CrankNicolson(thomas, A, B, U, n);
		}
		file << scientific << "  " << log10(errorMax(U, aU, n));

		//==================== KMB ====================
		U = initU(H[i], a, b, n);
		for (int k_ = 1; k_ < k_d; k_++)
		{
			KMB(U, lambda_d, n);
		}
		file << scientific << "  " << log10(errorMax(U, aU, n)) << endl;

		delete[] U;
		delete[] aU;
		delete[] B;
		for (int i = 0; i < n; i++)
		{
			delete[] A[i];
		}
		delete[] A;
	}
	file.close();

	//===================================== POINT 2 =====================================
	h = 0.05;
	n = fabs(a - b) / h + 1;
	//for KMB
	dt_d = (lambda_d * h * h) / D;
	k_d = (tmax / dt_d) + 2;
	//for Crank-Nicolson
	dt_ind = (lambda_ind * h * h) / D;
	k_ind = (tmax / dt_ind) + 2;
	//U init
	double *U2 = new double[n];
	//analyU
	double *aU2 = new double[n];
	//init b
	double *B2 = new double[n];

	int tk1_ind = 240; //t = 0.6
	int tk2_ind = 480; //t = 1.2
	int tk3_ind = 720; //t = 1.8

	int tk1_d = 600;  //t = 0.6
	int tk2_d = 1200; //t = 1.2
	int tk3_d = 1800; //t = 1.8

	//==================== Matrix Crank-Nicolson ====================
	double **A2 = new double *[n];
	for (int i = 0; i < n; i++)
	{
		A2[i] = new double[n];
	}
	CrankNicolson_A(A2, lambda_ind, ALPHA, BETA, FI, PSI, n, h);

	//==================== Crank-Nicolson with SOR ====================
	U2 = initU(h, a, b, n);
	for (int k_ = 1; k_ < k_ind; k_++)
	{
		CrankNicolson_B(B2, THETA, GAMMA, U2, n, lambda_ind);
		CrankNicolson(SOR, A2, B2, U2, n);
		if (tk1_ind == k_) //t == 0.6
		{
			analyU(aU2, h, a, b, n, 0.6);
			file.open("CN_SOR_t1.txt", ios::out);
			file << "x\tCN_SOR\tANALITYCZNE" << endl;
			for (int i = 0; i < n; i++)
			{
				file << scientific << a + i * h << " " << U2[i] << " " << aU2[i] << endl;
			}
			file.close();
		}
		if (tk2_ind == k_) //t == 1.2
		{
			analyU(aU2, h, a, b, n, 1.2);
			file.open("CN_SOR_t2.txt", ios::out);
			file << "x\tCN_SOR\tANALITYCZNE" << endl;
			for (int i = 0; i < n; i++)
			{
				file << scientific << a + i * h << " " << U2[i] << " " << aU2[i] << endl;
			}
			file.close();
		}
		if (tk3_ind == k_) //t == 1.8
		{
			analyU(aU2, h, a, b, n, 1.8);
			file.open("CN_SOR_t3.txt", ios::out);
			file << "x\tCN_SOR\tANALITYCZNE" << endl;
			for (int i = 0; i < n; i++)
			{
				file << scientific << a + i * h << " " << U2[i] << " " << aU2[i] << endl;
			}
			file.close();
		}
	}

	//==================== Crank-Nicolson with Thomas ====================
	U2 = initU(h, a, b, n);
	for (int k_ = 1; k_ < k_ind; k_++)
	{
		CrankNicolson_B(B2, THETA, GAMMA, U2, n, lambda_ind);
		CrankNicolson(thomas, A2, B2, U2, n);
		if (tk1_ind == k_) //t == 0.6
		{
			analyU(aU2, h, a, b, n, 0.6);
			file.open("CN_THOMAS_t1.txt", ios::out);
			file << "x\tTHOMAS\tANALITYCZNE" << endl;
			for (int i = 0; i < n; i++)
			{
				file << scientific << a + i * h << " " << U2[i] << " " << aU2[i] << endl;
			}
			file.close();
		}
		if (tk2_ind == k_) //t == 1.2
		{
			analyU(aU2, h, a, b, n, 1.2);
			file.open("CN_THOMAS_t2.txt", ios::out);
			file << "x\tTHOMAS\tANALITYCZNE" << endl;
			for (int i = 0; i < n; i++)
			{
				file << scientific << a + i * h << " " << U2[i] << " " << aU2[i] << endl;
			}
			file.close();
		}
		if (tk3_ind == k_) //t == 1.8
		{
			analyU(aU2, h, a, b, n, 1.8);
			file.open("CN_THOMAS_t3.txt", ios::out);
			file << "x\tTHOMAS\tANALITYCZNE" << endl;
			for (int i = 0; i < n; i++)
			{
				file << scientific << a + i * h << " " << U2[i] << " " << aU2[i] << endl;
			}
			file.close();
		}
	}

	//==================== KMB ====================
	U2 = initU(h, a, b, n);
	for (int k_ = 1; k_ < k_d; k_++)
	{
		KMB(U2, lambda_d, n);
		if (tk1_d == k_) //t == 0.6
		{
			analyU(aU2, h, a, b, n, 0.6);
			file.open("KMB_t1.txt", ios::out);
			file << "x\tKMB\tANALITYCZNE" << endl;
			for (int i = 0; i < n; i++)
			{
				file << scientific << a + i * h << " " << U2[i] << " " << aU2[i] << endl;
			}
			file.close();
		}
		if (tk2_d == k_) //t == 1.2
		{
			analyU(aU2, h, a, b, n, 1.2);
			file.open("KMB_t2.txt", ios::out);
			file << "x\tKMB\tANALITYCZNE" << endl;
			for (int i = 0; i < n; i++)
			{
				file << scientific << a + i * h << " " << U2[i] << " " << aU2[i] << endl;
			}
			file.close();
		}
		if (tk3_d == k_) //t == 1.8
		{
			analyU(aU2, h, a, b, n, 1.8);
			file.open("KMB_t3.txt", ios::out);
			file << "x\tKMB\tANALITYCZNE" << endl;
			for (int i = 0; i < n; i++)
			{
				file << scientific << a + i * h << " " << U2[i] << " " << aU2[i] << endl;
			}
			file.close();
		}
	}

	delete[] U2;
	delete[] aU2;
	delete[] B2;
	for (int i = 0; i < n; i++)
	{
		delete[] A2[i];
	}
	delete[] A2;

	//===================================== POINT 3 =====================================
	h = 0.05;
	n = fabs(a - b) / h + 1;
	//for KMB
	dt_d = (lambda_d * h * h) / D;
	k_d = (tmax / dt_d) + 2;
	//for Crank-Nicolson
	dt_ind = (lambda_ind * h * h) / D;
	k_ind = (tmax / dt_ind) + 2;
	//U init
	double *U3 = new double[n];
	//analyU
	double *aU3 = new double[n];
	//init b
	double *B3 = new double[n];

	//==================== Matrix Crank-Nicolson ====================
	double **A3 = new double *[n];
	for (int i = 0; i < n; i++)
	{
		A3[i] = new double[n];
	}
	CrankNicolson_A(A3, lambda_ind, ALPHA, BETA, FI, PSI, n, h);

	//==================== Crank-Nicolson with SOR ====================
	file.open("maxErrorCN_SOR.txt", ios::out);
	file << "maxerr\tt" << endl;
	U3 = initU(h, a, b, n);
	analyU(aU3, h, a, b, n, 0.0);
	file << errorMax(U3, aU3, n) << "  " << dt_ind * 0.0 << endl;
	for (int k_ = 1; k_ < k_ind; k_++)
	{
		analyU(aU3, h, a, b, n, dt_ind * k_);
		CrankNicolson_B(B3, THETA, GAMMA, U3, n, lambda_ind);
		CrankNicolson(SOR, A3, B3, U3, n);
		file << errorMax(U3, aU3, n) << "  " << dt_ind * k_ << endl;
	}
	file.close();

	//==================== Crank-Nicolson with Thomas ====================
	file.open("maxErrorCN_THOMAS.txt", ios::out);
	file << "maxerr\tt" << endl;
	U3 = initU(h, a, b, n);
	analyU(aU3, h, a, b, n, 0.0);
	file << errorMax(U3, aU3, n) << "  " << dt_ind * 0.0 << endl;
	for (int k_ = 1; k_ < k_ind; k_++)
	{
		analyU(aU3, h, a, b, n, dt_ind * k_);
		CrankNicolson_B(B3, THETA, GAMMA, U3, n, lambda_ind);
		CrankNicolson(thomas, A3, B3, U3, n);
		file << errorMax(U3, aU3, n) << "  " << dt_ind * k_ << endl;
	}
	file.close();

	//==================== KMB ====================
	file.open("maxErrorKMB.txt", ios::out);
	file << "maxerr\tt" << endl;
	U3 = initU(h, a, b, n);
	analyU(aU3, h, a, b, n, 0.0);
	file << errorMax(U3, aU3, n) << "  " << dt_d * 0.0 << endl;
	for (int k_ = 1; k_ < k_d; k_++)
	{
		analyU(aU3, h, a, b, n, dt_d * k_);
		KMB(U3, lambda_d, n);
		file << errorMax(U3, aU3, n) << "  " << dt_d * k_ << endl;
	}
	file.close();

	delete[] U3;
	delete[] aU3;
	delete[] B3;
	for (int i = 0; i < n; i++)
	{
		delete[] A3[i];
	}
	delete[] A3;

	return 0;
}