#include<math.h>
#include<stdio.h>
#include<stdlib.h>

#include <iostream>
#include <fstream>
#include <random>
#include <vector>
#include <math.h>

#include <string.h>

#include <numeric>
#include <complex>
#include <cmath>
#include <iomanip>



// -- Random number -----------
#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define NTAB 32
#define NDIV (1+(IM-1)/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)
#define PI acos(-1.0)
// ----------------------------


#define NR_END 1
#define FREE_ARG char*
//#define NE 3


using namespace std;









class PCOD{

	
public:
	struct particle_{
	double r_x;
	double r_y;
	double theta;
	double w;
	int id;

	complex<double> ck;
	complex<double> ck_current;

	complex<double> *ck_tau;  // Array with the past values of c_k
	double *theta_tau; // Array with the past values of theta

	// Used by the RK4, because we need to retrieve the values at: t - tau + h/2. 
	// Such values are not known, them one needs to interpolate known past values.
	complex<double> interp_ck_tau ; // Interpolated past value of c_k(t - tau + h/2)
	double interp_theta_tau; // Interpolated past value of theta(t - tau + h/2)
	};


	float ran1(long *idum_);

	void init(int N_, double M_, double omega0_, double h_);
	void step_forward();
	void destroy();

	int get_ticks();
	double get_t();


	particle_ *particles;
	double *d_theta;

	double K = 0.1;

	

private:
	


	// Numerical integration
	void nrerror(char error_text[]);
	void derivs(double y[],double df[], particle_ *particles, double interp_t, bool store);
	double *dvector(long nl,long nh);
	void free_dvector(double *v, long nl, long nh);

	
	// Interpolation
	void polint(double xa[], double ya[], int n, double x, double *y, double *dy);
	void retrieve_past_values_via_interpolation(particle_ particles[], double interp_t, int num_points);


	// Model
	void projectionMatrix();
	complex<double> scalar_product(double *v1, complex<double> *v2);
	double *matrix_line(double *matrix, int l);
	void rotation_center(double *state, particle_ *particles);
	void rotation_center_delay(double *state, particle_ *particles, bool store);
	double distance_particles(particle_ *particles, int id1, int id2);

	double uk_circular_symmetric_paley_all_to_all(particle_ *particles, double state[], double theta, int id);
	double uk_circular_symmetric_paley_all_to_all_delay(particle_ *particles, double theta, int id);


	

	// Random number generation Seed
	long idum{-986967};


	// --- Runge-Kutta 4th order --------
	double h=0.1;
	double tf = 5000.0;

	double t;
	int it;

	
	double *x, *a, *b, *c,*df,*y;
	double Pi_2 = 2.0*M_PI; 
	//-----------------------------------


	// ----- Model parameters ------------
	int N = 12;
	double M = 3;

	int NE = 36;

	double K_m = 0.18;
	double K_M = -0.02;

	
	double K0 = 0.1;
	double omega0 = 0.05; //0.05;

	double *P;
	complex<double> *cc;

	// delay
	bool delay;
	double tau;
	int delay_array_size = 100;
	int delay_array_size_minus_1 = 99;

	double h_over_2;
	double *t_array;
	int idx_tau; // Index for the arrays with the past values of tau and ck

	// -----------------------------------

};



