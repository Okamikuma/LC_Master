#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "hsph2d.h"

// Simple program aimed at checking the proper behavior of "hsph2d.h" functions.


// I think I can give default values in the header, it's just like object-oriented programming

int main(){
    
//  Global variables values are set here
    Nsph = 8;
    eta = 0.644;
    Ncol = 0;
    
//    sigma = 2.0f*sqrt(eta/((double)Nsph*PI));
    
    initsys();

	return 0;

}
