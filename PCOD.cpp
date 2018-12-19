#include "PCOD.h"



// ------------- Random number generator -----------------------------------

// Uniform distribution
/*Minimal random number generator of Park and Miller with Bays-Durham shuffle and added safeguards. Returns a uniform random deviate between 0.0 and 1.0 (exclusive of the endpoint values). Call with idum a negative integer to initialize; thereafter, do not alter idum
between successive deviates in a sequence. RNMX should approximate the largest floating value that is less than 1.*/
float PCOD::ran1(long *idum_){


	int j;
	long k;
	static long iy=0;
	static long iv[NTAB];
	float temp;

	if(*idum_<=0 || !iy){
		if(-(*idum_)<1) 
			*idum_=1;
		else 
			*idum_ = -(*idum_);
		for(j=NTAB+7;j>=0;j--){
			k=(*idum_)/IQ;
			*idum_=IA*(*idum_-k*IQ)-IR*k;
			if(*idum_<0) 
				*idum_ +=IM;
			if(j<NTAB) 
				iv[j]=*idum_;
		}
		iy=iv[0];
	}
	k=(*idum_)/IQ;
	*idum_=IA*(*idum_-k*IQ)-IR*k;
	if(*idum_<0) 
		*idum_ += IM;
	j=iy/NDIV;
	iy=iv[j];
	iv[j]=*idum_;
	if((temp=AM*iy)>RNMX) 
		return RNMX;
	else 
		return temp;
}



// ------------------- MODEL ---------------------------


void PCOD::projectionMatrix(){

	P = new double[N*N];

	for(int l=0; l<N; ++l){
		//std::cout << std::endl;
		for(int c=0; c<N; ++c){
			//std::cout << l << "  " << c << std::endl;
			if(l == c){
				*(P + l*N + c) = (double(N)-1.0) / double(N);
			}else{
				*(P + l*N + c) = -1.0 / double(N);
			}

			//std::cout << *(P + l*N + c) << "  ";
		}		
	}
}



complex<double> PCOD::scalar_product(double *v1, complex<double> *v2){
	complex<double> dot = 0.0;

	for(int i=0; i<N; ++i){
		dot += *(v1+i) * *(v2+i);
	}
	
	return dot;
}

double *PCOD::matrix_line(double *matrix, int l){
	double *Pk = new double[N];


	for(int i=0; i<N; ++i)
		*(Pk+i) = *(matrix + l*N + i);

	return Pk;
}



void PCOD::rotation_center(double *state, particle_ *particles){

	double r_x, r_y, theta;
	
	int index_c = 0;
	for(int i=0; i<NE; i+=3){

		theta = state[i+1];
		r_x = state[i+2];
		r_y = state[i+3];
		
		complex<double> rk(r_x,r_y);
		complex<double> vel(cos(theta),sin(theta));
		complex<double> ck( std::real(rk)- (std::imag(vel)/ particles[index_c].w), std::imag(rk) + ( std::real(vel) / particles[index_c].w));

		cc[index_c] = ck;

		++index_c;
	}
}


void PCOD::rotation_center_delay(double *state, particle_ *particles, bool store){

	double r_x, r_y, theta;
	int i;
	int index_c = 0;

	for(i=1; i<=NE; i+=3){
		// cc is the center of rotation in the past: (t-tau) or (t-tau+h/2) or (t-tau+h)...		
		cc[index_c] = particles[index_c].interp_ck_tau;

		theta = state[i];
		r_x = state[i+1];
		r_y = state[i+2];

		complex<double> rk(r_x,r_y);
		complex<double> vel(cos(theta),sin(theta));
		complex<double> ck( std::real(rk)- (std::imag(vel)/ particles[index_c].w), std::imag(rk) + ( std::real(vel) / particles[index_c].w));

		// Center at time t (current)
		particles[index_c].ck_current = ck;

		if(store){
			// Store the last (truly) numerically integrated ck in the array 
			particles[index_c].ck_tau[idx_tau] = particles[index_c].ck;
			cc[index_c] = particles[index_c].ck_tau[idx_tau];
			particles[index_c].interp_ck_tau = particles[index_c].ck_tau[idx_tau];

			particles[index_c].ck = ck;
		}
	
		++index_c;
	}
}


double PCOD::distance_particles(particle_ *particles, int id1, int id2){
	double d = 0.0;
	
	d = sqrt( (particles[id1].r_x - particles[id2].r_x) * (particles[id1].r_x - particles[id2].r_x) + 
	    (particles[id1].r_y - particles[id2].r_y) * (particles[id1].r_y - particles[id2].r_y) );

	return d;
}


double PCOD::uk_circular_symmetric_paley_all_to_all(particle_ *particles, double state[], double theta, int id){

	double *Pk = matrix_line(P, id);

	complex<double> eitk(cos(theta), sin(theta));
	complex<double> Pk_c = scalar_product(Pk, cc);

	double potential_derivative = 0.0;
	double Km = 0.0;
	double theta_j=0.0;
	
	int m,j;	

	for(m=1; m<M+1; ++m){
		if(m == M)
			Km = K_M;
		else
			Km = K_m;
		
		for(j=1; j<=NE; j+=3){
			theta_j = state[j];
			potential_derivative += (Km/double(m)) * sin(double(m) * ( theta_j - theta ));	
		}		
	}


	potential_derivative = potential_derivative / double(N);

	double ukk = particles[id].w * (1.0 + K * std::real(std::conj(eitk) * Pk_c)) - potential_derivative;

	delete Pk;

	return ukk;
}




double PCOD::uk_circular_symmetric_paley_all_to_all_delay(particle_ *particles, double theta, int id){

	double *Pk = matrix_line(P, id);

	complex<double> eitk(cos(theta), sin(theta));
	complex<double> Pk_c = scalar_product(Pk, cc);

	double potential_derivative = 0.0;
	double Km = 0.0;
	double theta_j=0.0;

	int m, j;
	
	for(m=1; m<M+1; ++m){
		if(m == M)
			Km = K_M;
		else
			Km = K_m;
	
		for(j=0; j<N; ++j){

			// Every j except j=id, i.e., no self-delay.
			if(j != id){
				theta_j = particles[j].interp_theta_tau;
				potential_derivative += (Km/double(m)) * sin(double(m) * ( theta_j - theta ));	
			}
		}		
	}

	potential_derivative = potential_derivative / double(N);

	double ukk = particles[id].w * (1.0 + K * std::real(std::conj(eitk) * Pk_c)) - potential_derivative;

	delete Pk;

	return ukk;
}




// ------------------------ END MODEL --------------------------------------







// --------- Numerical Integration -------------------------------
void PCOD::nrerror(char error_text[]){
	fprintf(stderr,"Numerical Recipes run-time error...\n");
	fprintf(stderr,"%s\n",error_text);
	fprintf(stderr,"...now exiting to system...\n");
	exit(1);
}

void PCOD::derivs(double y[],double df[], particle_ *particles, double interp_t, bool store){

	int index=0, i;

	if(delay){
		// Claculate the rotation center of the particles
		rotation_center_delay(y, particles, store);

		for(i=1; i<=NE; i+=3){
			// The vector cc of centers contain delayed centers. 
			// We replace the c_k by its instantaneous center, because we are not considering self-delay.
			cc[index] = particles[index].ck_current;

			df[i]=uk_circular_symmetric_paley_all_to_all_delay(particles, y[i], index);
			df[i+1]=cos(y[i]);
			df[i+2]=sin(y[i]);

			// Replacing the instantaneous c_k by its delayed version for future calculations.
			cc[index] = particles[index].interp_ck_tau;

			++index;
		}
	}else{
		rotation_center(y, particles);

		for(i=1; i<=NE; i+=3){
			df[i]=uk_circular_symmetric_paley_all_to_all(particles, y, y[i], index);
			df[i+1]=cos(y[i]);
			df[i+2]=sin(y[i]);

			++index;
		}
	}
}

double *PCOD::dvector(long nl,long nh){
	double *v;

	v=(double *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(double)));
	return v-nl+NR_END;
}


void PCOD::free_dvector(double *v, long nl, long nh){
	free((FREE_ARG) (v+nl-NR_END));
}
// --------------------------------------------------------------------



void PCOD::polint(double xa[], double ya[], int n, double x, double *y, double *dy){
// Adapted from Numerical Recipes: http://www.it.uom.gr/teaching/linearalgebra/NumericalRecipiesInC/c3-1.pdf
// Given arrays xa[1..n] and ya[1..n], and given a value x, this routine returns a value y
// and an  error estimate dy. If P(x) is the polynomial of degree N − 1 such that P(xai)= yai; i =1 ;:::; n,
// then  the  returned  value y=P(x).

	int i,m,ns=1;
	double den,dif,dift,ho,hp,w;

	double *c,*d;
	dif=fabs(x-xa[1]);
	c=dvector(1,n);
	d=dvector(1,n);
	for (i=1;i<=n;i++) {
		// Here we  find the index ns of the closest table entry,
		if ( (dift=fabs(x-xa[i])) < dif) {
			ns=i;
			dif=dift;
		}
		c[i]=ya[i];
		// and  initialize  the tableau of c's and d's.
		d[i]=ya[i];
	}
	*y=ya[ns--];

	// This  is  the  initial  approximation  to y.
	for (m=1;m<n;m++) {
		//For  each  column  of the  tableau,
		for (i=1;i<=n-m;i++) {
			// we  loop  over  the  current c's  and d's  and  update them.
			ho=xa[i]-x;
			hp=xa[i+m]-x;
			w=c[i+1]-d[i];
			if ( (den=ho-hp) == 0.0) nrerror("Error in routine polint");
			//This error can occur only if  two input xa's are (to within  roundo  ) identical.
			den=w/den;
			d[i]=hp*den;
			//Here  the c's and d's are updated.
			c[i]=ho*den;
		}
		*y += (*dy=(2*ns < (n-m) ? c[ns+1] : d[ns--]));
		// After  each  column  in  the  tableau  is  completed,  we  decide  which  correction, c or d,
		// we want  to  add  to our  accumulating  value  of y,  i.e., which  path to take  through the
		// tableau|forking up or down.  We do this in such a way as to take the most "straight
		// line"  route through the  tableau to its  apex, updating ns accordingly to keep track of
		// where we are.  This route keeps the partial approximations centered (insofar as possible)
		// on the target x. Thelast dy added  is  thus the  error  indication.
	}
	free_dvector(d,1,n);
	free_dvector(c,1,n);
}





void PCOD::retrieve_past_values_via_interpolation(particle_ particles[], double interp_t, int num_points){
	// We have theta(t) and want to calculate theta(t+h)
	// We have stored the values [ theta(t-tau-h), theta(t-tau), theta(t-tau+h), theta(t-tau+2h), .., theta(t-h)]
	// The Runge-Kutta 4th algorithm for the delayed case needs the value theta(t-tau+h/2). 
	// Such values are not stored, then one needs to interpolate neighboring known values of theta.
	// Here, we are using the 4 neghbors [ theta(t-tau-h), theta(t-tau), theta(t-tau+h), theta(t-tau+2h) ]

	int array_lim = num_points + 1;

	// Arrays for the interpolation algorithm
	double *time = dvector(1,array_lim);
	double *cx_tau = dvector(1,array_lim);
	double *cy_tau = dvector(1,array_lim);
	double *theta_tau = dvector(1,array_lim);

	//double time[num_points];
	//double cx_tau[num_points];
	//double cy_tau[num_points];
	//double theta_tau[num_points];
	double error;

	double cx_tau_;
	double cy_tau_;

	int index=0;
	int index_ini;
	int i, j;


	index_ini = idx_tau;
	if(index_ini > delay_array_size_minus_1)
			index_ini = 0;	

	for(i=0; i<N; ++i){

		// Initializing the index for the delayed values array
		index = index_ini;

		// Retrieve the 4 known values around t-tau
		for(j=1; j<array_lim; ++j){

			time[j] = t_array[index];
			cx_tau[j] = std::real(particles[i].ck_tau[index]);
			cy_tau[j] = std::imag(particles[i].ck_tau[index]);
			theta_tau[j] = particles[i].theta_tau[index];

			++index;
			if(index > delay_array_size_minus_1)
				index = 0;		
		}

		

		// Calculate the interpolation of such values to estimate theta_k and ck at t-tau+h/2 or t-tau+h

		// -------------------- theta(interp_t) = .. ----------------------------------------------
		// The specified time is given by the argument interp_t. 
		polint(time, theta_tau, num_points, interp_t, &particles[i].interp_theta_tau, &error);

		// If the error is greater than 10^3, we use linear interpolation
		if(error > fabs(0.001))
			particles[i].interp_theta_tau = (theta_tau[1] + theta_tau[2]) / 2.0;

		particles[i].interp_theta_tau = atan2(sin(particles[i].interp_theta_tau), cos(particles[i].interp_theta_tau));
		// ---------------------------------------------------------------------------------------


		// ---------------------- ck(interp_t) ----------------------------------------------------
		polint(time, cx_tau, num_points, interp_t, &cx_tau_, &error);
		// If the error is greater than 10^3, we use linear interpolation
		if(error > fabs(0.001))
			cx_tau_ = (cx_tau[1] + cx_tau[2]) / 2.0;

		polint(time, cy_tau, num_points, interp_t, &cy_tau_, &error);
		// If the error is greater than 10^3, we use linear interpolation
		if(error > fabs(0.001))
			cy_tau_ = (cy_tau[1] + cy_tau[2]) / 2.0;

		complex<double> ck(cx_tau_, cy_tau_);
		particles[i].interp_ck_tau = ck;
		// ---------------------------------------------------------------------------------------



		/*
		// Linear interpolation -----------------
		// One uses it to faster simulations, by the cost of an increasing error

		particles[i].interp_theta_tau = (theta_tau[1] + theta_tau[2]) / 2.0;

		cx_tau_ = (cx_tau[1] + cx_tau[2]) / 2.0;
		cy_tau_ = (cy_tau[1] + cy_tau[2]) / 2.0;
		complex<double> ck(cx_tau_, cy_tau_);
		particles[i].interp_ck_tau = ck;
		// ---------------------------------------
		*/


		// Check if this is a NaN
		if(particles[i].interp_theta_tau != particles[i].interp_theta_tau || particles[i].interp_ck_tau != particles[i].interp_ck_tau){
			cout << "Either theta_tau or ck_tau is NaN" << endl;
			for(j=0; j<delay_array_size; ++j){
				cout << "j:" << j << "   theta_tau=" << particles[0].theta_tau[j] << "    ck=" <<  particles[0].ck_tau[j] << endl;
			}		
			exit(1);

		}
	}

	// Free memory
	free_dvector(time,1,array_lim);
	free_dvector(cx_tau,1,array_lim);
	free_dvector(cy_tau,1,array_lim);
	free_dvector(theta_tau,1,array_lim);
}







int PCOD::get_ticks(){
	return it;
}

double PCOD::get_t(){
	return t;
}



void PCOD::init(int N_, double M_, double omega0_, double tau_){
	int i,j;	

	N = N_;
	M = M_;
	omega0 = omega0_;
	tau = tau_;

	NE = N_ * 3;
	
	K_m = 0.18;
	K_M = -0.02;

	K = 0.3;

	particles = new particle_[N]; // (particle_ *)malloc((size_t) ((N_)*sizeof(particle_)));

	y=dvector(1,NE+2);
	df=dvector(1,NE+2);

	x=dvector(1,NE+1);
	a=dvector(1,NE+2);
	b=dvector(1,NE+2);
	c=dvector(1,NE+2);

	for(i=1; i<=NE; ++i)
		x[i] = 1;


	// -------------- Delay -------------------------------------------
	if(tau > 0)
		delay = true;
	else
		delay = false;

	// It is divided by delay_array_size because this is the size of the delay array.
	// The integration step should be a divisor of tau
	h = tau_ / double(delay_array_size);
	if(h == 0.0)
		h = 0.1;

	h_over_2 = h/2.0;

	// Creating an array with all time values in [t-tau,t)
	t_array = new double[delay_array_size];
	for(i=0; i<delay_array_size; ++i)
		t_array[i] = -tau + i*h;

	// ----------------------------------------------------------------

	idx_tau = 0;
	t=0.0;
	it=0;

	double box_side = 15.0; //4.0*sqrt(double(N)*M_PI); //1.0/(3.0*omega0);

	cc = new complex<double>[N];
	projectionMatrix();
	d_theta = new double[N];

	
	// Initial conditions
	for(i=0; i<N; ++i){
		particles[i].id = i;
		particles[i].theta = ran1(&idum) * 2.0 * M_PI;
		particles[i].w = omega0;

		if(ran1(&idum) > 0.5)
			particles[i].r_x = ran1(&idum) * box_side;
		else
			particles[i].r_x = -ran1(&idum) * box_side;



		if(ran1(&idum) > 0.5)
			particles[i].r_y = ran1(&idum) * box_side;
		else
			particles[i].r_y = -ran1(&idum) * box_side;


		// Rotation center ck
		complex<double> rk(particles[i].r_x,particles[i].r_y);
		complex<double> vel(cos(particles[i].theta),sin(particles[i].theta));
		complex<double> ck( std::real(rk)- (std::imag(vel)/ particles[i].w), std::imag(rk) + ( std::real(vel) / particles[i].w));

		// Allocating memory
		particles[i].theta_tau = new double[delay_array_size];
		particles[i].ck_tau = new complex<double>[delay_array_size];

		// Past values of theta and ck. As we do not know, we choose the current values.
		for(j=0; j<delay_array_size; ++j){
			particles[i].theta_tau[j] = particles[i].theta; 
			particles[i].ck_tau[j] = ck;
		}
	
		particles[i].interp_theta_tau = particles[i].theta; // Interpolated theta at (t-tau+h/2) or (t-tau+h)
		particles[i].interp_ck_tau    = ck; // Interpolated center at (t-tau+h/2) or (t-tau+h)
		particles[i].ck = ck;  // Last calculated center at t multiple of h.
		particles[i].ck_current = ck; // Current center
		cc[i] = ck;

		d_theta[i] = omega0;		
	}
}



void PCOD::step_forward(){
	
	int j,i;
   
	double t_minus_tau = t - tau;
	
	int index = 0;
	for(i=1; i<=NE; i+=3){
		x[i]=particles[index].theta;
		x[i+1]=particles[index].r_x;	
		x[i+2]=particles[index].r_y;

		++index;
	}

	for(j=1;j<=NE;j++)
			y[j]=x[j];

	//cout << "DERIVS K1" << endl;
	// Depends on the states at (t-tau). We already have such values:
	for(i=0; i<N; ++i){
		particles[i].interp_theta_tau = particles[i].theta_tau[idx_tau];
		particles[i].interp_ck_tau = particles[i].ck_tau[idx_tau];
	}
	derivs(y,df, particles, t_minus_tau, true); // Delay at (t-tau)

	// k1 is a
	for(j=1;j<=NE;j++){
		a[j]=h*df[j];
		y[j]=x[j]+a[j]/2.0;
	}

	//cout << "DERIVS K2" << endl;
	// Needs interpolation
	// As the past values of theta and ck are not available, we interpolate.
	retrieve_past_values_via_interpolation(particles, t_minus_tau + h_over_2, 5);
	derivs(y,df, particles, t_minus_tau + h_over_2, false); // Delay at (t-tau+h/2)

	// k2 is b
	for(j=1;j<=NE;j++){
		b[j]=h*df[j];
		y[j]=x[j]+b[j]/2.0;
	}

	//cout << "DERIVS K3" << endl;
	// Needs interpolation at the same point as before, then we maintain the values we got
	derivs(y,df, particles, t_minus_tau + h_over_2, false); // Delay at (t-tau+h/2)

	// k3 is c
	for(j=1;j<=NE;j++){
		c[j]=h*df[j];
		y[j]=x[j]+c[j];
	}

	//cout << "DERIVS K4" << endl;
	// Noes not need interpolation, because we already have the states at t_=t-tau+h:
	index=idx_tau+1;
	if(index > delay_array_size_minus_1)
		index = 0;	
	for(j=0; j<N; ++j){
		particles[j].interp_theta_tau = particles[j].theta_tau[index];
		particles[j].interp_ck_tau = particles[j].ck_tau[index];
	}
	derivs(y,df, particles, t_minus_tau+h, false); // Delay at (t-tau+h)

	// k4 is h*df[j]
	for(j=1;j<=NE;j++)
		x[j]=x[j]+(a[j]+h*df[j])/6.0+(b[j]+c[j])/3.0;

	index = 0;
	for(i=1; i<=NE; i+=3){

		if(delay){
			// Store the past value of theta(t-tau) before replacing the variable by theta(t)
			particles[index].theta_tau[idx_tau] = particles[index].theta;
		}

		particles[index].theta = x[i];
		particles[index].r_x = x[i+1];
		particles[index].r_y = x[i+2];

		// Updating the output for the ALF
		d_theta[index] = (  (a[i]+h*df[i])/6.0+(b[i]+c[i])/3.0  ) / h;

		if(particles[index].theta > Pi_2)
			particles[index].theta -= Pi_2;
		else if(particles[index].theta < 0.0)
			particles[index].theta += Pi_2;

		++index;
	}

	if(delay){
		t_array[idx_tau] += tau;// t_minus_tau;	

		++idx_tau;
		if(idx_tau > delay_array_size_minus_1) 
			idx_tau = 0;

	}		
	
	t=t+h;
	it++;
}



void PCOD::destroy(){

	// Free memmory
	delete P;
	delete cc;
	delete t_array;
	delete d_theta;
	free_dvector(y,1,NE+2);
	free_dvector(df,1,NE+2); 

	free_dvector(x,1,NE+1);
	free_dvector(a,1,NE+2);
	free_dvector(b,1,NE+2);
	free_dvector(c,1,NE+2);
}





