#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "random.h"



int main(){
	double f = -756.323980992309834673;
	double *v;
	int i;

	v= (double *)malloc(10*sizeof(double));
	rlxd_init(1,1);
	
	ranlxd(v,10);
	for(i=0;i<10;i++)
	{
		v[i]*=f;
	        printf("\n %lf %lf \n", floor(v[i]),v[i]);

	}
	
	printf("\n %lf \n", floor(f));
	return 0;

}
