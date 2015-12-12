#include <iostream>
#include <cmath>

using namespace std;
// global variables
const double TOL = 1E-5;
const double mu = 0.012277471;

// Define functions
void k_func(double Y4[], double k[][7], double dt, double a[][7], double c[]);
double k_0_calc(double k_vec, double x1, double x2,double x3,double x4);
double k_1_calc(double k_vec, double x1, double x2,double x3,double x4);
double k_2_calc(double k_vec, double x1, double x2,double x3,double x4);
double k_3_calc(double k_vec, double x1, double x2,double x3,double x4);
void rk_4(double Y4[], double k[][7], double dt);
void rk_5(double Y5[], double k[][7], double dt);
double new_step(double Y4[], double Y5[], double dt);

// Main function
int main(){
	double dt = 1E-5;
	// Time variable
	double T = 0;
	double T_end = 12.0;

	// For RK4 -> Y4 = [ x, y, x', y']
	double* Y4 = new double [4];

	// For RK5 -> Y5 = [ x, y, x', y']
	double* Y5 = new double [4];

    double k[4][7];

	// Because I alreay know the size of Matrix A
	// I use static allocation as multidimensional Array.
	// Only the lower diagonal matrix will be != 0.
	// Yeah, I know we should use a long array instead a "matrix" in C++
	// But in this case we already know the dimension.
    double a[7][7] =
    {
    { 0,            0,              0,              0,          0,              0,       0 },
    { 1.0/5,        0,              0,              0,          0,              0,       0 },
    { 3.0/40,       9.0/40,         0,              0,          0,              0,       0 },
    { 44.0/45,      -56.0/15,       32.0/9,         0,          0,              0,       0 },
    { 19372.0/6561, -25360.0/2187,  64448.0/6561,   -212.0/729, 0,              0,       0 },
    { 9017.0/3168,  -355.0/33,      46732.0/5247,   49.0/176,   -5103.0/18656,  0,       0 },
    { 35.0/384,     0,              500.0/1113,     125.0/192,  -2187.0/6784,   11.0/84, 0 }
    };

	// The c-vector is for both rk methods the same
	// Because of that, I define it one time only in the main
	// function.
	double c[7] = {0, 1.0/5,  3.0/10, 4.0/5,  8.0/9,  1,  1 };

	// Initial values at T=0
	Y4[0] = 0.994;
    Y4[1] = 0;
	Y4[2] = 0;
    Y4[3] = -2.00158510637908;

	Y5[0] = 0.994;
	Y5[1] = 0;
    Y5[2] = 0;
	Y5[3] = -2.00158510637908;
	cout << T << "\t" << Y4[0] << "\t" << Y4[1] << endl;

	while(T < T_end){
        // Calculate k vectors
        k_func(Y4, k, dt, a, c);
        // Use specific weights for RK4 method
        rk_4(Y4, k, dt);
        // Use specific weights for RK5 method
        rk_5(Y5, k, dt);
        // Give out new y_n+1
        cout << T << "\t" << Y4[0] << "\t" << Y4[1] << endl;
        //cout << T << "\t" << Y5[0] << "\t" << Y5[1] << endl;
        // One step more
        T += dt;
        dt = new_step(Y4, Y5, dt);
        Y5[0] = Y4[0];
        Y5[1] = Y4[1];
        Y5[2] = Y4[2];
        Y5[3] = Y4[3];
	}
	return 0;
}

void k_func(double Y4[], double k[][7], double dt, double a[][7], double c[]){
    double tmp_0, tmp_1, tmp_2, tmp_3;

	for(int i=0; i<7; i++){
    tmp_0 = 0;
    tmp_1 = 0;
    tmp_2 = 0;
    tmp_3 = 0;
    for(int j=0; j<i; j++){
        tmp_0 += a[i][j]*k[i-1][0];
        tmp_1 += a[i][j]*k[i-1][1];
        tmp_2 += a[i][j]*k[i-1][2];
        tmp_3 += a[i][j]*k[i-1][3];
        //cout << "a(" << i << ";" << j << ")*k[" << i-1 << "] ";
    }
    //cout << endl;
    k[i][0] = k_0_calc(k[i][0], Y4[0]+dt*tmp_0, Y4[1]+dt*tmp_1, Y4[2]+dt*tmp_2, Y4[3]+dt*tmp_3);
    k[i][1] = k_1_calc(k[i][1], Y4[0]+dt*tmp_0, Y4[1]+dt*tmp_1, Y4[2]+dt*tmp_2, Y4[3]+dt*tmp_3);
    k[i][2] = k_2_calc(k[i][2], Y4[0]+dt*tmp_0, Y4[1]+dt*tmp_1, Y4[2]+dt*tmp_2, Y4[3]+dt*tmp_3);
    k[i][3] = k_3_calc(k[i][3], Y4[0]+dt*tmp_0, Y4[1]+dt*tmp_1, Y4[2]+dt*tmp_2, Y4[3]+dt*tmp_3);
    cout << k[i][0] << "\t" << k[i][1] << "\t" << k[i][2] << "\t" << k[i][3] << endl;
	}
}

double k_0_calc(double k_vec, double x1, double x2,double x3,double x4){
	k_vec = x3;
	return k_vec;
}

double k_1_calc(double k_vec, double x1, double x2,double x3,double x4){
	k_vec = x4;
	return k_vec;
}

double k_2_calc(double k_vec, double x1, double x2,double x3,double x4){
    double r = sqrt(pow(x1+mu,2)+pow(x2,2));
    double s = sqrt(pow(x1-1+mu,2)+pow(x2,2));
	k_vec = x1 + 2.0*x4 - (1.0-mu)*(x1+mu)/pow(r,3)  - mu*(x1-1.0+mu)/pow(s,3);
	return k_vec;
}

double k_3_calc(double k_vec, double x1, double x2,double x3,double x4){
    double r = sqrt(pow(x1+mu,2)+pow(x2,2));
    double s = sqrt(pow(x1-1+mu,2)+pow(x2,2));
	k_vec = x2 - 2.0*x3 - (1.0-mu)*x2/pow(r,3)       - mu*x2/pow(s,3);
	return k_vec;
}

void rk_4(double Y4[], double k[][7], double dt){
	// Define specific constant weights for rk_4
	double* b = new double [7];

	b[0] = 5179.0/57600;
	b[1] = 0;
	b[2] = 7571.0/16695;
	b[3] = 393.0/640;
	b[4] = -92097.0/339200;
	b[5] = 187.0/2100;
	b[6] = 1.0/40;
    double tmp_0 = 0;
    double tmp_1 = 0;
    double tmp_2 = 0;
    double tmp_3 = 0;
	for(int i=0 ; i<7 ; i++)
	{
		tmp_0 += b[i]*k[0][i];
		tmp_1 += b[i]*k[1][i];
		tmp_2 += b[i]*k[2][i];
		tmp_3 += b[i]*k[3][i];
	}
    Y4[0] += dt*tmp_0;
	Y4[1] += dt*tmp_1;
    Y4[2] += dt*tmp_2;
	Y4[3] += dt*tmp_3;
	//cout << Y4[0] << "\t" << Y4[1] << "\t" << Y4[2] << "\t" << Y4[3] << endl;
}

void rk_5(double Y5[], double k[][7], double dt){
	// Define specific constant weights for rk_5
	double* b = new double [7];

	b[0] = 35.0/384;
	b[1] = 0;
	b[2] = 500.0/1113;
	b[3] = 125.0/192;
	b[4] = -2187.0/6784;
	b[5] = 11.0/84;
	b[6] = 0;
    double tmp_0 = 0;
    double tmp_1 = 0;
    double tmp_2 = 0;
    double tmp_3 = 0;

	for(int i=0 ; i<7 ; i++)
	{
		tmp_0 += b[i]*k[0][i];
		tmp_1 += b[i]*k[1][i];
		tmp_2 += b[i]*k[2][i];
		tmp_3 += b[i]*k[3][i];
	}
    Y5[0] += dt*tmp_0;
	Y5[1] += dt*tmp_1;
    Y5[2] += dt*tmp_2;
	Y5[3] += dt*tmp_3;
}

double new_step(double Y4[], double Y5[], double dt){
	// Security Factor q:=0.5
	double q = 0.5;
	// Order of the RK method
	int p = 4;
	// Calculate H (local error)
	double H = max(
        abs(Y4[0]-Y5[0]),
        abs(Y4[1]-Y5[1]));
    H = max(
        H,
        abs(Y4[2]-Y5[2]));
    H = max(
        H,
        abs(Y4[3]-Y5[3]));

	// Calculate new dt
	//dt = dt*q*pow( (TOL/H) , (1.0/(p+1)) );
	dt = 1E-5;
	return dt;
}
