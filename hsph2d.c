#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include "random.h"

#define PI 3.141592653589793

/*Disposizione BCC per le coordinate delle sfere*/

void BCCdisp(double *cx, double *cy, int n, double d)
{
	int i1, i2, i3, lato;
	double passo;
	
	i1 = 0;
	i2 = 0;
	i3 = 0;
	lato = 0;
	passo = 0.0;
	
	/*Le sfere vengono divise in due gruppi che andranno a formare il reticolo e il reticolo traslato diagonale. Ognuno dei due reticoli viene preso completo per questo si richiede un numero di sfere che sia il doppio di un quadrato perfetto.
        Definizione del numero di sfere per lato (poi tutto il reticolo verrà raddoppiato traslato in diagonale)*/
    
	lato = sqrt(n/2);
	
	/*Definizione del passo del reticolo bcc*/
	passo = (1.0 - d)/( (double)(lato) - 0.5);
	
	/*Disposizione dei reticoli*/
	for( i1 = 0; i1 < lato; i1++ )
	{
		for( i2 = 0; i2 < lato; i2++ )
		{
			*(cx + i3) = i2*passo + d/2.0;
			
			*(cy + i3) = i1*passo + d/2.0;
			
			*(cx + i3 + lato*lato) = *(cx + i3) + passo/2.0; 
			
			*(cy + i3 + lato*lato) = *(cy + i3) + passo/2.0;
			
			i3++;
		}
	}
}

/*Funzione che calcola il prodotto scalare tra i due vettori (a1, a2) e (b1, b2) (in due dimensioni)*/
double prodscal( double a1, double a2, double b1, double b2 )
{
	double ps;
	
	ps = 0.0;
	
	ps = (a1*b1) + (a2*b2);
	
	return ps;
}

/*Funzione che calcola la distanza tra i centri di due sfere*/
double distanza( double *cx, double *cy, int n1, int n2 )
{
	double dist;
	
	dist = 0.0;
	
	dist = sqrt( prodscal( *(cx + n1) - *(cx + n2), *(cy + n1) - *(cy + n2), *(cx + n1) - *(cx + n2), *(cy + n1) - *(cy + n2) ) );
	
	return dist;
}

/*Funzione che calcola il tempo di collisione tra due sfere date le loro posizioni e velocità (e tenendo conto delle copie)*/
double tempocollisione( double *cx, double *cy, double *velx, double *vely, int ns, double d, int n1, int n2 )
{
	double t, r12[2], v12[2], discr, tempocopia[9], coordcopia[2], r12v12, mqv, mqr;
	int i1, i2, i3;
	
	t = 0.0;
	discr = 0.0;
	r12v12 = 0.0;
	mqv = 0.0;
	mqr = 0.0;
	i1 = 0;
	i2 = 0;
	i3 = 0;
	
	for( i1 = 0; i1 < 2; i1++ )
	{
		r12[i1] = 0.0;
		v12[i1] = 0.0;
		coordcopia[i1] = 0.0;
	}
	
	for( i1 = 0; i1 < 9; i1++ )
	{
		tempocopia[i1] = 0.0;
	}
	
	/*Si setta in partenza t = -1.0 (nessuna collisione)*/
	t = -1.0;
	
	/*Con gli indici i1 e i2 si scorre sulle copie della particella n2, si calcolano i tempi di collisione e, se c'è un tempo di collisione diverso da -1.0 (infinito), allora quello è il tempo di collisione tra la sfera n1 e la sfera n2*/
	for( i1 = -1; i1 < 2; i1++ )
	{
		for( i2 = -1; i2 < 2; i2++ )
		{
			coordcopia[0] = *(cx + n2) + (double)(i2);
			coordcopia[1] = *(cy + n2) + (double)(i1);
			
			r12[0] = *(cx + n1) - coordcopia[0];
			r12[1] = *(cy + n1) - coordcopia[1];
			
			v12[0] = *(velx + n1) - *(velx + n2);
			v12[1] = *(vely + n1) - *(vely + n2);
			
			r12v12 = prodscal(r12[0], r12[1], v12[0], v12[1]);
			
			mqv = prodscal(v12[0], v12[1], v12[0], v12[1]);
			
			mqr = prodscal(r12[0], r12[1], r12[0], r12[1]);
			
			discr = r12v12*r12v12 - mqv*( mqr - d*d );
			
			if( ( r12v12 >= 0.0 )||( discr < 0.0 ) )
			{
				tempocopia[i3] = -1.0;
			}
			
			else
			{
				tempocopia[i3] = - r12v12 - sqrt( discr );
				
				tempocopia[i3] = tempocopia[i3]/mqv;
				
				if(tempocopia[i3] < 0.0)
				{
					printf("\nAttezione: tempo negativo. prodscal = %lf, discr = %lf, r12v12^2 = %lf, sqrt(mqr) = %lf, d = %lf\n", r12v12, discr, r12v12*r12v12, sqrt(mqr), d);
				}
			}
			
			i3 = i3 + 1;
		
			coordcopia[0] = 0.0;
			coordcopia[1] = 0.0;
			discr = 0.0;
			r12[0] = 0.0;
			r12[1] = 0.0;
			v12[0] = 0.0;
			v12[1] = 0.0;
			mqv = 0.0;
			mqr = 0.0;
			r12v12 = 0.0;
		}
	}
	
	/*Nel vettore tempocopia[i] ci sono i tempi di collisione di n1 con tutte le copie di n2. Il tempo di collisione tra n1 ed n2 sarà il minore tempo positivo nel vettore (oppure -1.0 se non ve ne è)*/
	for( i1 = 0; i1 < 9; i1++ )
	{
		if( tempocopia[i1] != -1.0 )
		{
			if( t == -1.0 )
			{
				t = tempocopia[i1];
			}
			
			else
			{
				if( tempocopia[i1] < t )
				{
					t = tempocopia[i1];
				}
			}
		}
	}
	
	return t;
}

/*Funzione che cerca e restituisce il tempo minimo (tra quelli positivi) nella matrice dei tempi*/
double tempominimo( double **matr, int ns, int *p1, int *p2 )
{
	double tmin;
	int i1, i2;
	
	tmin = -1.0;
	i1 = 0;
	i2 = 0;
	
	for( i1 = 1; i1 < ns; i1++ )
	{
		for( i2 = 0; i2 < i1; i2++ )
		{
			if( *(*(matr + i1) + i2) != -1.0 )
			{
				if( tmin == -1.0 )
				{
					tmin = *(*(matr + i1) + i2);
					
					*p1 = i1;
					*p2 = i2;
				}
				
				else
				{
					if( *(*(matr + i1) + i2) < tmin )
					{
						tmin = *(*(matr + i1) + i2);
						
						*p1 = i1;
						*p2 = i2;
					}
				}
			}
		}
	}
	
	return tmin;
}

/*Funzione che restituisce il modulo di un numero*/
double modulo( double p )
{
	double modp;
	
	modp = 0.0;
	
	if( p >= 0.0 )
	{
		modp = p;
	}
	
	else
	{
		modp = -p;
	}
	
	return modp;
}

/*Funzione che fa evolvere il sistema fino all'urto successivo, aggiorna le velocità delle particelle che collidono, aggiorna la matrice dei tempi e restituisce il tempo trascorso*/
double urto( double *coordinate_x, double *coordinate_y, double *velocita_x, double *velocita_y, double **matrice, int numero_sfere, double d, double *dvelx, double *dvely )
{
	double tm, congiungente[2], vel_relativa[2], ausiliario, delta[2];
	int indice1, indice2, part1, part2;
	
	tm = 0.0;
	ausiliario = 0.0;
	indice1 = 0;
	indice2 = 0;
	
	for( indice1 = 0; indice1 < 2; indice1++ )
	{
		congiungente[indice1] = 0.0;
		vel_relativa[indice1] = 0.0;
		delta[indice1] = 0.0;
	}
	
	/*Calcolo il tempo minimo nella matrice*/
	tm = tempominimo( matrice, numero_sfere, &part1, &part2 );	
	
	/*Faccio evolvere il sistema con moto rettilineo uniforme*/
	for( indice1 = 0; indice1 < numero_sfere; indice1++ )
	{
		*(coordinate_x + indice1) = *(coordinate_x + indice1) + *(velocita_x + indice1)*tm;
		*(coordinate_y + indice1) = *(coordinate_y + indice1) + *(velocita_y + indice1)*tm;
	}
	
	/*Impongo le condizioni di bordo periodiche*/
	for( indice1 = 0; indice1 < numero_sfere; indice1++ )
	{
		while( *(coordinate_x + indice1) > 1.0 )
		{
			*(coordinate_x + indice1) = *(coordinate_x + indice1) - 1.0;
		}
		
		while( *(coordinate_x + indice1) < 0.0 )
		{
			*(coordinate_x + indice1) = *(coordinate_x + indice1) + 1.0;
		}
		
		while( *(coordinate_y + indice1) > 1.0 )
		{
			*(coordinate_y + indice1) = *(coordinate_y + indice1) - 1.0;
		}
		
		while( *(coordinate_y + indice1) < 0.0 )
		{
			*(coordinate_y + indice1) = *(coordinate_y + indice1) + 1.0;
		}
	}
	
	/*Sottraggo il tempo dell'urto alle entrate della matrice diverse da -1.0*/
	for( indice1 = 1; indice1 < numero_sfere; indice1++ )
	{
		for( indice2 = 0; indice2 < indice1; indice2++ )
		{
			if( *(*(matrice + indice1) + indice2) != -1.0 )
			{
				*(*(matrice + indice1) + indice2) = *(*(matrice + indice1) + indice2) - tm;
			}
		}
	}
	
	/*Tengo conto delle copie*/
	if( modulo( *(coordinate_x + part1) - *(coordinate_x + part2) ) > 0.5 )
	{
		if( *(coordinate_x + part1) - *(coordinate_x + part2) > 0.0 )
		{
			delta[0] = 1.0;
		}
		
		if( *(coordinate_x + part1) - *(coordinate_x + part2) < 0.0 )
		{
			delta[0] = -1.0;
		}
	}
		
	if( modulo( *(coordinate_y + part1) - *(coordinate_y + part2) ) > 0.5 )
	{
		if( *(coordinate_y + part1) - *(coordinate_y + part2) > 0.0 )
		{
			delta[1] = 1.0;
		}
		
		if( *(coordinate_y + part1) - *(coordinate_y + part2) < 0.0 )
		{
			delta[1] = -1.0;
		}
	}
	
	/*Calcolo le componenti del versore congiungente e del vettore velocità relativa*/
	congiungente[0] = *(coordinate_x + part1) - (*(coordinate_x + part2) + delta[0]);
	congiungente[1] = *(coordinate_y + part1) - (*(coordinate_y + part2) + delta[1]);
	
	ausiliario = prodscal( congiungente[0], congiungente[1], congiungente[0], congiungente[1] );
	ausiliario = sqrt(ausiliario);
	
	congiungente[0] = congiungente[0]/ausiliario;
	congiungente[1] = congiungente[1]/ausiliario;
	
	ausiliario = 0.0;
	
	vel_relativa[0] = *(velocita_x + part1) - *(velocita_x + part2);
	vel_relativa[1] = *(velocita_y + part1) - *(velocita_y + part2);
	
	/*Metto l'opposto delle componenti della velocità relativa prima dell'urto in dvelx e dvely, in modo tale che sommando quelle dopo l'urto si ottenga la variazione della velocità relativa*/
	*dvelx = -vel_relativa[0];
	*dvely = -vel_relativa[1];
	
	/*Aggiorno le componenti delle velocità*/
	ausiliario = prodscal( vel_relativa[0], vel_relativa[1], congiungente[0], congiungente[1] );
	
	*(velocita_x + part1) = *(velocita_x + part1) - ausiliario*congiungente[0];
	*(velocita_y + part1) = *(velocita_y + part1) - ausiliario*congiungente[1];
	*(velocita_x + part2) = *(velocita_x + part2) + ausiliario*congiungente[0];
	*(velocita_y + part2) = *(velocita_y + part2) + ausiliario*congiungente[1];
	
	/*Calcolo le componenti della velocità relativa DOPO l'urto*/
	vel_relativa[0] = *(velocita_x + part1) - *(velocita_x + part2);
	vel_relativa[1] = *(velocita_y + part1) - *(velocita_y + part2);
	
	/*Calcolo la variazione delle componenti della velocità relativa prima e dopo l'urto (finale - iniziale)*/
	*dvelx = *dvelx + vel_relativa[0];
	*dvely = *dvely + vel_relativa[1];
	
	ausiliario = prodscal( vel_relativa[0], vel_relativa[1], congiungente[0], congiungente[1] );
	
	/*Aggiorno le entrate della matrice nelle righe/colonne corrispondenti alle particelle che hanno colliso*/
	for( indice1 = 0; indice1 < numero_sfere; indice1++ )
	{
		*(*(matrice + part1) + indice1) = tempocollisione( coordinate_x, coordinate_y, velocita_x, velocita_y, numero_sfere, d, part1, indice1 );
	}
	
	for( indice1 = 0; indice1 < numero_sfere; indice1++ )
	{
		*(*(matrice + indice1) + part1) = tempocollisione( coordinate_x, coordinate_y, velocita_x, velocita_y, numero_sfere, d, indice1, part1 );
	}
	
	for( indice1 = 0; indice1 < numero_sfere; indice1++ )
	{
		*(*(matrice + part2) + indice1) = tempocollisione( coordinate_x, coordinate_y, velocita_x, velocita_y, numero_sfere, d, part2, indice1 );
	}
	
	for( indice1 = 0; indice1 < numero_sfere; indice1++ )
	{
		*(*(matrice + indice1) + part2) = tempocollisione( coordinate_x, coordinate_y, velocita_x, velocita_y, numero_sfere, d, indice1, part2 );
	}
		
	return tm;
}

/*Funzione che restituisce la parte intera (come int) di un numero double*/
int parteintera( double g )
{
	int pint;
	
	pint = 0;
	
	pint = (int)(g);
	
	return pint;
}

/*Velocities histogram function*/
void hist(double *sample, int N, int Nb, int Nc, FILE *F_fr)
{
    double smp_max, amp_b, aux, cen_b, norm;
    int i;
    double *frqs;
    //    inizializzazione delle variabili dichiarate
    smp_max=0.0;
    amp_b=0.0;
    aux=0.0;
    i = 0;
    frqs = NULL;
    cen_b = 0.0;
    norm = 0.0;
    //    Trovo la velocità massima in modulo e poi definisco l'ampiezza del bin
    for( i = 0; i< N*Nc; i++ )
    {
        if( modulo(*(sample + i)) > modulo(smp_max) )
            smp_max = *(sample + i);
    }
    //  smp_max is varied so that maximum values will be slightly smaller and will fall into the bin without problems
    smp_max = smp_max*(1.00000001)
    /*L'ampiezza di ogni bin è 2*maxv/Nbin*/
    amp_b = (2.0*modulo(smp_max))/((double)(Nb));
    
    
    /*Alloco memoria per i vettori in cui ciascuna posizione corrisponde a un bin e il valore contenuto è la frequenza relativa*/
    frqs = (double *) malloc(Nb*sizeof(double));
    
    /*Scorro il vettore dei dati e per ciascun dato aumento di 1 la componente del vettore delle frequenze corrispondente al bin in cui il dato capita*/
    for( i = 0; i < Nc*N; i++ )
        //  To the smp_max is added a small number in order to be bigger than the maximum.
        //  otherwise this code needs to be uncommented.Beware! This may affect heavily the distribution of bins.
        //        if( *(sample + i) == smp_max )
        //        {
        //            if( *(sample + i) > 0.0 )
        //            {
        //                *(frqs + Nb - 1) = *(frqs + Nb - 1) + 1.0;
        //            }
        //
        //            if( *(sample + i) < 0.0 )
        //            {
        //                *(frqs + 0) = *(frqs + 0) + 1.0;
        //            }
        //        }
        
    /*Per ogni valore della velocità si aumenta di 1 la frequenza del bin corrispondente. I bin devono essere pari. In tal caso una velocità v positiva capita nel bin numero (Nb/2 + parteintera(v/ampiezza_bin)), mentre se v < 0 capiterà nel bin numero ((Nb/2 - 1+ parteintera(v/ampiezza_bin)). La parte intera è con il suo segno*/
        else if( *(sample + i) >= 0.0 )
        {
            *(frqs + parteintera(*(sample + i)/amp_b) + Nb/2) = *(frqs + parteintera(*(sample + i)/amp_b) + Nb/2) + 1.0;
        }
    
        else if( *(sample + i) < 0.0 )
        {
            *(frqs + parteintera(*(sample + i)/amp_b) + Nb/2 - 1) = *(frqs + parteintera(*(sample + i)/amp_b) + Nb/2 - 1) + 1.0;
        }
    centrobin = - modulo(smp_max) + (amp_b/2.0);
    
    /*Al fine di ottenere il corretto istogramma (pdf con area 1), ogni frequenza relativa deve essere moltiplicata per una normalizzazione*/
    norm = 1.0/amp_b;
    
    for( i = 0; i < Nb; i++ )
    {
        /*Divido le frequenze per il numero totale di dati, in modo da avere le frequenze relative*/
        *(freqs + i) = ((*(freqs + i))/((double)(ns*nc)))*norm;
        
        /*Scrivo centro del bin e frequenza relativa del bin nei file*/
        fprintf( F_fr, "%lf %lf\n", centrobin[0], *(freqs + i) );
        centrobin = centrobin + amp_b;
        //          Ottengo alla fine due colonne nella prima c'è la posizione del bin e nella seconda la frequenza.
    }
    
    free( freqs );
}



int main(){
	int i, j, nsfere, nurti, nbin, k;
	double diametro, frimp, tempo, vmedia[2], encin, temperatura, press, auxilium, intervallofrimp;
	double *x, *y, *vx, *vy, *dativx, *dativy, **datitempicollisioni;
	double **matricetempi;
	FILE *filefreqvx, *filefreqvy;
	
	/*Inizializzo le variabili*/
	i = 0;
	j = 0;
	k = 0;
	nsfere = 0;
	nurti = 0;
	nbin = 0;
	diametro = 0.0;
	frimp = 0.0;
	tempo = 0.0;
	encin = 0.0;
	temperatura = 0.0;
	press = 0.0;
	auxilium = 0.0;
	intervallofrimp = 0.0;
	x = NULL;
	y = NULL;
	vx = NULL;
	vy = NULL;
	dativx = NULL;
	dativy = NULL;
	datitempicollisioni = NULL;
	matricetempi = NULL;
	
	for( i = 0; i < 2; i++ )
	{
		vmedia[i] = 0.0;
	}
	
	/*Titolo*/
	printf("\n\n__________SFERE RIGIDE IN 2 DIMENSIONI__________\n\n");
	
	/*Richiedo il numero di sfere*/
	printf("Numero sfere (usare 2, 8, 18, 32, 50, 72, 98, 128, 162, 200, ...): ");
	scanf("%d", &nsfere);
	
	/*Richiedo la frazione di impacchettamento*/
	printf("\n\nFrazione di impacchettamento: ");
	scanf("%lf", &frimp);
	
	/*Calcolo il diametro*/
	diametro = 2.0*sqrt(frimp/((double)(nsfere)*PI));
	printf("\n\nDiametro = %lf\n", diametro);
	
	/*Alloco memoria per posizioni e velocità*/
	x = (double *) malloc( nsfere*(sizeof(double)) );
	y = (double *) malloc( nsfere*(sizeof(double)) );
	vx = (double *) malloc( nsfere*(sizeof(double)) );
	vy = (double *) malloc( nsfere*(sizeof(double)) );
	
	/*Inizializzo le posizioni delle sfere*/
	for( i = 0; i < nsfere; i++ )
	{
		*(x + i) = 0.0;
		*(y + i) = 0.0;
		*(vx + i) = 0.0;
		*(vy + i) = 0.0;
	}
	
	/*Creo la disposizione bcc*/
	BCCdisp( x, y, nsfere, diametro );
	
	/*Controllo che le sfere non entrino una nell'altra a causa della frazione di impacchettamento troppo grande*/
	for( i = 0; i < nsfere -1; i++ )
	{
		for( j = i + 1; j < nsfere; j++ )
		{		
			if( distanza(x, y, i, j) <= diametro )
			{
				printf("\n\nAttenzione: la frazione di impacchettamento è troppo alta e le sfere %d e %d entrano una nell'altra\n\n", i, j);
				
				printf("\ndistanza = %lf\n", distanza(x, y, i, j));
				
				return 0;
			}
		}
	}

	/*Inizializzazione generatore numeri casuali*/
	rlxd_init(1, 1);

	/*Generazione casuale delle velocità*/
	ranlxd( vx, nsfere );
	ranlxd( vy, nsfere );

	/*I numeri sono estratti in [0, 1]: raddoppio e shifto di -1 per avere velocità in [-1, 1].*/
	for( i = 0; i < nsfere; i++ )
	{
		*(vx + i) = 2.0*(*(vx + i)) - 1.0;

		*(vy + i) = 2.0*(*(vy + i)) - 1.0;
	}
	
	/*Calcolo le componenti della velocità media*/
	for( i = 0; i < nsfere; i++ )
	{
		vmedia[0] = vmedia[0] + *(vx + i);
		vmedia[1] = vmedia[1] + *(vy + i);
	}
	
	vmedia[0] = vmedia[0]/(double)(nsfere);
	vmedia[1] = vmedia[1]/(double)(nsfere);
	
	/*Shifto le velocità per avere momento totale nullo*/
	for( i = 0; i < nsfere; i++ )
	{
		*(vx + i) = *(vx + i) - vmedia[0];
		*(vy + i) = *(vy + i) - vmedia[1];
	}
	
	/*Alloco memoria per la matrice dei tempi*/
	matricetempi = (double **) malloc(nsfere*sizeof(double *));
	
	for( i = 0; i < nsfere; i++ )
	{
		*(matricetempi + i) = (double *) malloc(nsfere*sizeof(double));
	}
	
	/*Inizializzo le entrate della matrice*/
	for( i = 0; i < nsfere; i++ )
	{
		for( j = 0; j < nsfere; j++ )
		{
			*(*(matricetempi + i) + j) = 0.0;
		}
	}
	
	/*Riempio metà matrice dei tempi con i tempi di collisione di tutte le sfere, l'altra metà è simmetrica, e non sarà mai letta*/
	for( i = 0; i < nsfere; i++ )
	{
		for( j = 0; j <= i; j++ )
		{
			*(*(matricetempi + i) + j) = tempocollisione( x, y, vx, vy, nsfere, diametro, i, j );
		}
	}
	
	printf("\nTempo iniziale = %lf\n", tempo);
	
	printf("\n\nNumero urti per la termalizzazione: ");
	scanf("%d", &nurti);
	printf("\n\n");
	
	/*Si realizzano nurti urti*/
	for( i = 0; i < nurti; i++ )
	{	
		tempo = tempo + urto( x, y, vx, vy, matricetempi, nsfere, diametro, &auxilium, &auxilium );
	
		printf("Tempo = %lf\n", tempo);
	}
	
	/*Calcolo dell'energia cinetica totale*/
	for( i = 0; i < nsfere; i++ )
	{
		encin = encin + (*(vx + i))*(*(vx + i)) + (*(vy + i))*(*(vy + i));
	}
	
	encin = encin/2.0;
	
	/*Calcolo della temperatura dal teorema di equipartizione dell'energia*/
	temperatura = encin/(double)(nsfere);
	
	printf("\nEnergia cinetica totale = %lf\nTemperatura*kb = %lf\n", encin, temperatura );
	
	/*Acquisisco il numero di campionamenti delle velocità per l'istogramma. Il numero di dati sarà nsfere*nurti*/
	printf("\n\nNumero campionamenti per l'istogramma delle velocità: ");
	scanf("%d", &nurti);
	
	/*Acquisisco il numero di bin (deve essere pari)*/
	printf("\nNumero di bin (pari): ");
	scanf("%d", &nbin);
	printf("\n");
	
	/*Alloco memoria per i vettori che conterranno tutti i valori campionati delle velocità*/
	dativx = (double *) malloc(nsfere*nurti*sizeof(double));
	dativy = (double *) malloc(nsfere*nurti*sizeof(double));
	
	/*Si realizza un campionamento di velocità ogni nsfere urti*/
	for( i = 0; i < nurti; i++ )
	{
		for( j = 0; j < nsfere; j++ )
		{
			tempo = tempo + urto( x, y, vx, vy, matricetempi, nsfere, diametro, &auxilium, &auxilium );
		}
		
		printf("Tempo = %lf\n", tempo);
		
		for( j = 0; j < nsfere; j++ )
		{
			*(dativx + i*nsfere + j) = *(vx + j);
			*(dativy + i*nsfere + j) = *(vy + j);
		}
	}
	
	filefreqvx = fopen("fv_x.txt", "w");
	filefreqvy = fopen("fv_y.txt", "w");
    hist(dativx,nsfere,nbin,nurti,filefreqvx);
    hist(dativy,nsfere,nbin,nurti,filefreqvy);
		
	free( x );
	free( y );
	free( vx );
	free( vy );
	free( dativx );
	free( dativy );
	
	for( i = 0; i < nsfere; i++ )
	{
		free( *(matricetempi + i) );
	}
	
	free( matricetempi );
	
	fclose( filefreqvx );
	fclose( filefreqvy );

	return 0;
}
