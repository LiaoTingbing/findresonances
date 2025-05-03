// findresonances.cpp: 定义应用程序的入口点。
//

#include "../include/main.h"


using namespace std;

int main()
{
	
	double f1 = 1.765;
	double alpha1 = 0.005;
	double ampl1 = 1.3;
	double phase1 = 0.4;
	double f2 = 2.345;
	double alpha2 = 0.012;
	double ampl2 = 0.45;
	double phase2 = 1.234;

	vec t_long = linspace(0, 20, 201);
	vec t_short = t_long.rows(0,19);
	cx_vec signal_long = ampl1 * exp(-IU*(2.0 * PI * f1 * t_long - phase1)) % exp(-alpha1 * t_long) + ampl2 * exp(-IU*(2 * PI * f2 * t_long - phase2)) % exp(-alpha2 * t_long);
	//cx_vec signal_short = ampl1 * exp(-IU*(2.0 * PI * f1 * t_short - phase1)) % exp(-alpha1 * t_short) + ampl2 * cos(2 * PI * f2 * t_short - phase2) % exp(-alpha2 * t_short);
	//signal_short.print();
	findresonances(signal_long, t_long, vec{ 1.5,10 },4,1e-7 );
	return 0;
}
