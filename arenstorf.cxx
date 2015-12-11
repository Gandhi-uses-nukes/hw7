#include <iostream>
#include <cmath>

using namespace std;
// global variables
const double TOL = 1E-5;
const double mu = 0.012277471;

// Define functions
void k_func(double Y4_x[], double Y4_y[], double kx[], double ky[], double dt, double a[][7], double c[]);
double k_x_calc(double Y4_x[], double Y4_y[], double r, double s, double factor_x, double factor_y);
double k_y_calc(double Y4_x[], double Y4_y[], double r, double s, double factor_x, double factor_y);
void rk_4(double Y4_x[], double Y4_y[], double kx[], double ky[], double dt);
void rk_5(double Y5_x[], double Y5_y[], double kx[], double ky[], double dt);
double new_step(double Y4_x[], double Y4_y[], double Y5_x[], double Y5_y[], double dt);

// Main function
int main(){
	double dt = 10E-5;
	// Time variable
	double T = 0;
	double T_end = 18.0;

	// For RK4 -> Y4 = [ y(t_n),  y'(t_n), y''(t_n) ]
	double* Y4_x = new double [3];
	double* Y4_y = new double [3];
	// For RK5 -> Y5 = [ x, y, x', y']
	double* Y5_x = new double [4];
	// Y5_y = [ y(t_n),  y'(t_n), y''(t_n) ]
	double* Y5_y = new double [3];

    // Define the k-vector for x and y
    double* kx = new double [7];
    double* ky = new double [7];

	// define constans for the different RK methods
	double* c = new double [7];

	// Because I alreay know the size of Matrix A
	// I use static allocation as multidimensional Array.
	// Only the lower diagonal matrix will be != 0.
	// Yeah, I know we should use a long array instead a "matrix" in C++
	// But in this case we already know the dimension.
    double a[7][7];
    for(int i=0 ; i<7 ; i++){
        for(int j=0 ; j<7 ; j++)
            a[i][j] = 0;
    }
	a[1][0] = 1.0/5;

	a[2][0] = 3.0/40;
	a[2][1] = 9.0/40;

	a[3][0] = 44.0/45;
	a[3][1] = -56.0/15;
	a[3][2] = 32.0/9;

	a[4][0] = 19372.0/6561;
	a[4][1] = -25360.0/2187;
	a[4][2] = 64448.0/6561;
	a[4][3] = -212.0/729;

	a[5][0] = 9017.0/3168;
	a[5][1] = -355.0/33;
	a[5][2] = 46732.0/5247;
	a[5][3] = 49.0/176;
	a[5][4] = -5103.0/18656;

	a[6][0] = 35.0/384;
	a[6][1] = 0;
	a[6][2] = 500.0/1113;
	a[6][3] = 125.0/192;
	a[6][4] = -2187.0/6784;
	a[6][5] = 11.0/84;

	// The c-vector is for both rk methods the same
	// Because of that, I define it one time only in the main
	// function.
	c[0] = 0;
	c[1] = 1.0/5;
	c[2] = 3.0/10;
	c[3] = 4.0/5;
	c[4] = 8.0/9;
	c[5] = 1;
	c[6] = 1;

	// Initial values at T=0
	Y4_x[0] = 0.994;
    Y4_y[0] = 0;
	Y4_x[1] = 0;
    Y4_y[1] = -2.00158510637908;

	Y5_x[0] = 0.994;
	Y5_y[0] = 0;
    Y5_x[1] = 0;
	Y5_y[1] = -2.00158510637908;
	cout << T << "\t" << Y4_x[0] << "\t" << Y4_y[0] << endl;

	while(T < T_end){
        // Calculate k vector
        k_func(Y4_x, Y4_y, kx, ky, dt, a, c);
        // Use specific weights for RK4 method
        rk_4(Y4_x, Y4_y, kx, ky, dt);
        // Use specific weights for RK5 method
        rk_5(Y5_x, Y5_y, kx, ky, dt);
        // Give out new y_n+1
        cout << T << "\t" << Y4_x[0] << "\t" << Y4_y[0] << endl;
        // One step more
        T += dt;
        dt = new_step(Y4_x, Y4_y, Y5_x, Y5_y, dt);
	}
	return 0;
}

void k_func(double Y4_x[], double Y4_y[], double kx[], double ky[], double dt, double a[][7], double c[]){
    double factor_x = 0;
    double factor_y = 0;

	for(int i=0; i<7; i++){
		factor_x = 0;
		factor_y = 0;
		if (i!=0){
            for(int j=0; j<7; j++){
                factor_x += a[i][j]*dt*kx[i-1];
                factor_y += a[i][j]*dt*ky[i-1];
                //cout << "a(" << i << ";" << j << ") = " << a[i][j] << endl;
			}
		}

		double r = sqrt(pow((Y4_x[0]+factor_x)+mu,2) + pow((Y4_y[0]+factor_y),2));
		double s = sqrt(pow((Y4_x[0]+factor_x)-1+mu,2) + pow((Y4_y[0]+factor_y),2));
		kx[i] = k_x_calc(Y4_x, Y4_y, r, s, factor_x, factor_y);
		ky[i] = k_y_calc(Y4_x, Y4_y, r, s, factor_x, factor_y);
	}
}

double k_x_calc(double Y4_x[], double Y4_y[], double r, double s, double factor_x, double factor_y){
	return (Y4_x[0]+factor_x) + 2.0*(Y4_y[1]+factor_y) - (1.0-mu)*((Y4_x[0]+factor_x)+mu)/pow(r,3) - mu*((Y4_x[0]+factor_x)-1.0+mu)/pow(s,3);
}

double k_y_calc(double Y4_x[], double Y4_y[], double r, double s, double factor_x, double factor_y){
	return (Y4_y[0]+factor_y) - 2.0*(Y4_x[1]+factor_x) - (1.0-mu)*(Y4_y[0]+factor_y)/pow(r,3) - (Y4_x[0]+factor_x)/pow(s,3);
}

void rk_4(double Y4_x[], double Y4_y[], double kx[], double ky[], double dt){
	// Define specific constant weights for rk_4
	double* b = new double [7];

	b[0] = 5179.0/57600;
	b[1] = 0;
	b[2] = 7571.0/16695;
	b[3] = 393.0/640;
	b[4] = -92097.0/339200;
	b[5] = 187.0/2100;
	b[6] = 1/40;
    double tmp_x = 0;
    double tmp_y = 0;
	for(int i=0 ; i<7 ; i++)
	{
		tmp_x += b[i]*kx[i];
		tmp_y += b[i]*ky[i];
	}
    Y4_x[0] += dt*tmp_x;
	Y4_y[0] += dt*tmp_y;
}

void rk_5(double Y5_x[], double Y5_y[], double kx[], double ky[], double dt){
	// Define specific constant weights for rk_5
	double* b = new double [7];

	b[0] = 35.0/384;
	b[1] = 0;
	b[2] = 500.0/1113;
	b[3] = 125.0/192;
	b[4] = -2187.0/6784;
	b[5] = 11.0/84;
	b[6] = 0;
    double tmp_x = 0;
    double tmp_y = 0;
	for(int i=0 ; i<7 ; i++)
	{
		tmp_x += b[i]*kx[i];
		tmp_y += b[i]*ky[i];
	}
	Y5_x[0] += dt*tmp_x;
	Y5_y[0] += dt*tmp_y;
}

double new_step(double Y4_x[], double Y4_y[], double Y5_x[], double Y5_y[], double dt){
	// Security Factor q:=0.5
	double q = 0.5;
	// Order of the RK method
	int p = 4;
	// Calculate H (local error)
	double H = abs(Y4_x[0]-Y5_x[0]);
	// Calculate new dt
	dt = dt*q*pow( (TOL/H) , (1.0/(p+1)) );
	return dt;
}

