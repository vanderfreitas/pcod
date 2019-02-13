#include "PCOD.h"


using namespace std;




class particle_delay : public particle_{

public:
	complex<double> ck;
	complex<double> ck_current;

	complex<double> *ck_tau;  // Array with the past values of c_k
	double *theta_tau; // Array with the past values of theta

	// Used by the RK4, because we need to retrieve the values at: t - tau + h/2. 
	// Such values are not known, them one needs to interpolate known past values.
	complex<double> interp_ck_tau ; // Interpolated past value of c_k(t - tau + h/2)
	double interp_theta_tau; // Interpolated past value of theta(t - tau + h/2)
};




// The PCOD_delay class inherits from the PCOD class
class PCOD_delay : public PCOD{

	
public:

	/*
	// Overriding
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
	};*/

	PCOD_delay(int N_, double M_, double omega0_, double tau_);

	// Overriding method
	void step_forward();

	// Overriding method
	void destroy();

	particle_delay *particles;
	

private:

	// Overriding method
	void derivs(double y[],double df[], double interp_t, bool store);

	
	// Interpolation
	void polint(double xa[], double ya[], int n, double x, double *y, double *dy);
	void retrieve_past_values_via_interpolation(double interp_t, int num_points);


	// Model
	void rotation_center_delay(double *state, bool store);
	double uk_circular_symmetric_paley_all_to_all_delay(double theta, int id);
	
	// ----- Model attributes ------------
	bool delay;
	double tau;
	int delay_array_size = 100;
	int delay_array_size_minus_1 = 99;

	double h_over_2;
	double *t_array;
	int idx_tau; // Index for the arrays with the past values of tau and ck
	// -----------------------------------
};



