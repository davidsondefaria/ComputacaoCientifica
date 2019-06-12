//#include "ranf.h"
#include <time.h>
#include <math.h>

double dist(double med,double var){
	double u1=ranf();//(rand()%10000)/1e4;
	double ue= sqrt(-2.0 * log(u1)) * cos(2*M_PI * (rand()%10000)/1e4);
	//while(u1*var*sqrt(2*M_PI)>1)
	//	u1=ranf();
	//if(ranf()<0.5)
		//return (cos(2*M_PI*ranf())* (var*sqrt(-2*log(u1*var*sqrt(2*M_PI)) ) ) +med );
		return var*ue +med;
	//return (-1*var*sqrt(-2*log(u1*var*sqrt(2*M_PI)) ) +med );

}


