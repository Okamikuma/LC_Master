#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define PI 3.141592653589793

typedef sphere{
    // Position of each sphere in phase-space
    // int dim;
    double pos[2];
    //double vel[DIMENSION];
}sphere;

//  Initialization function
 // questa funzione deve ricevere tutto quello serve per definire velocita e posizione di una particella la farei che prende la lista di particelle e scrive direttamente le coordinate di queste particelle nello spazio delle fasi.

void initsys(int Nsph, sphere *sys, double sigma){
    
    // the first for runs through all the particles the second and the third are used to assign the values on the BCC
    // Nsph is constrained to be two times the square of an integer
    
    int id01, id02, id03, size;
    double step;
    double *x1;
    
    // allocation of memory
    x1 = (double *)malloc( 2*size*sizeof(double) );
    
    id01 = 0;
    id02 = 0;
    id03 = 0;
    size = 0;
    step = 0.0f;
    size = sqrt(Nsph/2);
    step = (1.0 - sigma)/((double)size - 0.5);
    
    // The coordinates are created
    for (id01=0; id01<size; id01++) {
        x1[id01] = id01*step + d/2.0;
    }
    
    // Coordinates are distributed to the spheres
    for(id01 = 0, id01 < Nsph, id01++)
    {
        if (id01 < Nsph/2) {
            for (id02 = 0; id02<size; id02++)
            {
                for (id03=0; id03<size; id03++)
                {
                    sys[id01].pos[1] = x1[id02];
                    sys[id01].pos[2] = x1[id03];
                }
            }
        }
        if (id02 >= Nsph/2) {
                for (id02 = 0; id02<size; id02++)
                {
                    for (id03=0; id03<size; id03++)
                    {
                        sys[id01].pos[1] = x1[id02]+step/2.0;
                        sys[id01].pos[2] = x1[id03]+step/2.0;
                    }
                }
        }
    }
 }

int main(){

    int Nsph01,id01;
    
    double eta;
    double sigma01;
    sphere *sys01;
    
    FILE *prova = fopen("prova.txt","w")
    double eta = 0.2;
    int Nsph01 = 200;
    
    sigma01 = 2.0*sqrt(eta/(double)Nsph01*PI);
    
    sys01 = (body*)malloc( Nsph01*sizeof(sphere) );
    
    initsys(Nsph01, sys01, sigma01);
    
    for(id01 = 0, id01 < Nsph01, id01++)
    {
        fprintf(prova,"%lf %lf\n",sys01[id01].pos[1],sys01[id01].pos[2]);
    }
    free( prova);
    return 0;
    
}