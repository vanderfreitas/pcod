#include<iostream>

#include "PCOD.h"



int main(int argc, char *argv[]){
	int N = 12;
	int M = 3;
	double omega0 = 0.05;
	double h = 0.1;

	int i;
	int kk=0;
	

	char name[50];
	name[0] = '\0';
	strcat(name, "out_trajectory.csv");
	ofstream arq_out(name);

	// Object that controls the model
	PCOD pcod_model;

	// Initiating the pcod model
	pcod_model.init(N, M, omega0, h);

	while(pcod_model.get_t() < 1500.0){
		// Integrate the model
		pcod_model.step_forward();

		// Exporting data
		if(pcod_model.get_ticks()%10 == 0){
			for(i=0; i<N; ++i){
				arq_out << pcod_model.get_ticks() << "\t" << i << "\t" << pcod_model.particles[i].r_x << "\t" << pcod_model.particles[i].r_y << "\t" << pcod_model.particles[i].theta << "\t" << M <<  "\n";
			}
		}
	}

	// Finalize the model
	pcod_model.destroy();	

	// Close output file
	arq_out.close();


	return 0;
}
