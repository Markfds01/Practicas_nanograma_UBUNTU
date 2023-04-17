#include "random_number.h"

unsigned int irr[256];
unsigned int ir1;
unsigned char ind_ran,ig1,ig2,ig3;
void ini_ran(int SEMILLA)/**Funcion para cambiar la semilla**/
{
    int INI,FACTOR,SUM,i;

    srand(SEMILLA);
    INI = SEMILLA;
    FACTOR = 67397;
    SUM = 7364893;

    for(i=0;i<256;i++)
    {
        INI = (INI*FACTOR+SUM);
        irr[i] = INI;
    }
    ind_ran = ig1 = ig2 = ig3 =0;
}
float Parisi_Rapuano(long *idum)
{
int j;
long k;
static long iy=0;
static long iv[NTAB];
float temp;
if (*idum <= 0 || !iy) { 				//Initialize.
    if (-(*idum) < 1) *idum=1; 			//Be sure to prevent idum = 0.
    else *idum = -(*idum);
    for (j=NTAB+7;j>=0;j--) { 			//Load the shuffle table (after 8 warm-ups).
        k=(*idum)/IQ;
        *idum=IA*(*idum-k*IQ)-IR*k;
        if (*idum < 0) *idum += IM;
        if (j < NTAB) iv[j] = *idum;
    }
    iy=iv[0];
}
k=(*idum)/IQ; 							//Start here when not initializing.
*idum=IA*(*idum-k*IQ)-IR*k; 			//Compute idum=(IA*idum) % IM without overflows by Schrage�s method.
if(*idum < 0) *idum += IM;
j=iy/NDIV; 								//Will be in the range 0..NTAB-1.
iy=iv[j]; 								//Output previously stored value and refill the  shuffle table.
iv[j] = *idum;
if ((temp=AM*iy) > RNMX) return RNMX; 	//Because users don�t expect endpoint values.
else return temp;
}

double rand_gaussiano(long* idum)
{
    double x = 0;
    double y = 0;
    double res = 0;
    x = Parisi_Rapuano(idum);
    y = Parisi_Rapuano(idum);
    res = sqrt(-2 * log(x)) * cos(2*PI*y);
    return res;

}