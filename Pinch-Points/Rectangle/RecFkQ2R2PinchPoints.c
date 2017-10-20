#include "stdlib.h"
#include "stdio.h"
#include "math.h"

#define SEED 7
#define Q 2
#define LENGTH 2001 // must be odd
#define WIDTH 1001 //must be odd
#define S 1048575//32767
#define OUTFILE "RecFkQ2R2PinchPointsSketch.ps"
#define OUTFILE12 "RecFkQ2R2PinchPoints12.txt"
#define OUTFILE14 "RecFkQ2R2PinchPoints14.txt"
#define OUTFILE34 "RecFkQ2R2PinchPoints34.txt"
#define OUTFILE23 "RecFkQ2R2PinchPoints23.txt"
#define OUTFILE1234 "RecFkQ2R2PinchPoints1234.txt"
#define RUNFILE "RunFile.txt"
#define PRINTFREQ 1000000
#define RUNSMAX 50000005

void    randinit(long seed);
void    find(long x, long y, long s, long t);
void	walk(void);

#define M 16383                     

#define GetFromStack0(X,Y) {X = stackx[gptr0 & S]; Y = stacky[gptr0 & S]; ++gptr0;}
#define PutOnStack0(X,Y)	{stackx[pptr0 & S]=X; stacky[pptr0 & S]=Y; ++pptr0; if  (pptr0 == gptr0) {printf("error40\n");break;}}

#define GetFromStack1(X,Y) {X = collide1x[gptr1 & S]; Y = collide1y[gptr1 & S]; ++gptr1;}
#define PutOnStack1(X,Y)	{collide1x[pptr1 & S]=X; collide1y[pptr1 & S]=Y; ++pptr1; if  (pptr1 == gptr1) {printf("error41\n");break;}}

#define GetFromStack2(X,Y) {X = collide2x[gptr2 & S]; Y = collide2y[gptr2 & S]; ++gptr2;}
#define PutOnStack2(X,Y)	{collide2x[pptr2 & S]=X; collide2y[pptr2 & S]=Y; ++pptr2; if  (pptr2 == gptr2) {printf("error42\n");break;}}

#define NewRandomInteger (++nd,ra[nd&M]=ra[(nd-471)&M]^ra[(nd-1586)&M]^ra[(nd-6988)&M]^ra[(nd-9689)&M])

long    ra[M+1], nd;
long	lat[LENGTH][WIDTH];
long	walk1[LENGTH][WIDTH];
long	walk2[LENGTH][WIDTH];
long	pinchpoint12[LENGTH][WIDTH];
long	pinchpoint14[LENGTH][WIDTH];
long	pinchpoint34[LENGTH][WIDTH];
long	pinchpoint23[LENGTH][WIDTH];
long	pinchpoint1234[LENGTH][WIDTH];
long 	collide1x[S+1], collide1y[S+1];
long 	collide2x[S+1], collide2y[S+1];
long 	stackx[S+1], stacky[S+1];
long    runs, prob, connectivity;
long	gptr0, pptr0, gptr1, pptr1, gptr2, pptr2, gptr3, pptr3, gptr4, pptr4;

FILE	*fpsketch;

int main(void) 
{
	
	long  x, y, xo, yo, k, m, centerx, centery;
    FILE  *fp12, *fp14, *fp34, *fp23, *fp1234, *fprunfile;
    
    randinit(SEED);
	
    prob = (int) (sqrt(Q)/(1+sqrt(Q)) * 2147483648.0 );
	
	centerx = (LENGTH-1)*1.0/2;
	centery = (WIDTH-1)*1.0/2+1;
	
	for(x = 1; x < LENGTH-1; ++x)
		for (y = 0; y < WIDTH; ++y)
		{
			walk1[x][y] = walk2[x][y] = 0; //set walk variables on medial lattice sites to 0
			if ((x == 1) || (x == LENGTH-2)) lat[x][y] = Q;	//wire left/right sides
			else if ((y == 0) || (y == WIDTH-1)) lat[x][y] = 0;	//wire dual top/bottom sides
			else lat[x][y] = Q + (int)((NewRandomInteger/2147483648.0)*Q);
			pinchpoint12[x][y] = pinchpoint14[x][y] = pinchpoint34[x][y] = pinchpoint23[x][y] = pinchpoint1234[x][y] = 0;
		}
	
	for (runs = 1; runs <= RUNSMAX; ++runs)
	{
		
		gptr1 = pptr1 = gptr2 = pptr2 = gptr3= pptr3 = gptr4 = pptr4 = 0;		
		
		for(x = 1; x < LENGTH-1; ++x)
			for (y = 1; y < WIDTH-1; ++y)
				if (lat[x][y] >= Q) lat[x][y] -= Q; //initialize regular lattice sites with new spins
		
		for (xo = 1; xo < LENGTH-1; xo = xo+2)
			for (yo = 1; yo < WIDTH-1; yo = yo+2)
				if ((m = lat[xo][yo]) < Q)
				{ 
					gptr0 = pptr0 = 0;
					PutOnStack0(xo,yo)
					k = Q + (int)((NewRandomInteger/2147483648.0)*Q);  
					do {
						GetFromStack0(x,y)
						if (lat[x][y] < Q) find(x, y, k, m);
					}while (gptr0 != pptr0);
				}
		
		walk();
		
		while (gptr1 != pptr1) {
			GetFromStack1(x,y)
			if (connectivity == 1) ++pinchpoint12[x][y];
			if (connectivity == 2) ++pinchpoint14[x][y];
		}
		
		while (gptr2 != pptr2) {
			GetFromStack2(x,y)
			if (connectivity == 1) ++pinchpoint34[x][y];
			if (connectivity == 2) ++pinchpoint23[x][y];
		}
		
		if (runs%PRINTFREQ == 0) {
			
			fprunfile = fopen(RUNFILE, "w");
			fprintf(fprunfile, "%15ld\n", runs);
			fclose(fprunfile);
			
			fp12 = fopen(OUTFILE12,"w");
			fprintf(fp12, "%15ld%15ld%15ld\n", centerx, centery, pinchpoint12[centerx][centery]);
			for (x = 1; x < LENGTH-1; ++x)
				for (y = 1; y < WIDTH-1; ++y)
					if ((x%2 == 1 && y%2 == 0) || (x%2 == 0 && y%2 == 1)) fprintf(fp12, "%15ld%15ld%15ld\n", x, y, pinchpoint12[x][y]);
			fclose(fp12);
			
			fp14 = fopen(OUTFILE14,"w");
			fprintf(fp12, "%15ld%15ld%15ld\n", centerx, centery, pinchpoint14[centerx][centery]);
			for (x = 1; x < LENGTH-1; ++x)
				for (y = 1; y < WIDTH-1; ++y)
					if ((x%2 == 1 && y%2 == 0) || (x%2 == 0 && y%2 == 1)) fprintf(fp14, "%15ld%15ld%15ld\n", x, y, pinchpoint14[x][y]);
			fclose(fp14);
			
			fp34 = fopen(OUTFILE34,"w");
			fprintf(fp12, "%15ld%15ld%15ld\n", centerx, centery, pinchpoint34[centerx][centery]);
			for (x = 1; x < LENGTH-1; ++x)
				for (y = 1; y < WIDTH-1; ++y)
					if ((x%2 == 1 && y%2 == 0) || (x%2 == 0 && y%2 == 1)) fprintf(fp34, "%15ld%15ld%15ld\n", x, y, pinchpoint34[x][y]);
			fclose(fp34);
			
			fp23 = fopen(OUTFILE23,"w");
			fprintf(fp12, "%15ld%15ld%15ld\n", centerx, centery, pinchpoint23[centerx][centery]);
			for (x = 1; x < LENGTH-1; ++x)
				for (y = 1; y < WIDTH-1; ++y)
					if ((x%2 == 1 && y%2 == 0) || (x%2 == 0 && y%2 == 1)) fprintf(fp23, "%15ld%15ld%15ld\n", x, y, pinchpoint23[x][y]);
			fclose(fp23);
			
			fp1234 = fopen(OUTFILE1234,"w");
			fprintf(fp12, "%15ld%15ld%15ld\n", centerx, centery, pinchpoint1234[centerx][centery]);
			for (x = 1; x < LENGTH-1; ++x)
				for (y = 1; y < WIDTH-1; ++y)
					if ((x%2 == 1 && y%2 == 0) || (x%2 == 0 && y%2 == 1)) fprintf(fp1234, "%15ld%15ld%15ld\n", x, y, pinchpoint1234[x][y]);
			fclose(fp1234);
			
		}
		
	}/* for runs */
	
}


void find(long x, long y, long s, long t)
{
    int		j, sideflag;
	long	dx[4] = {1, 0, -1,  0}, dy[4] = {0, 1, 0, -1};
	long	xp, yp;
	
	lat[x][y] = s; //update spin
	
    for (j = 0; j < 4; ++j)
    {  
        sideflag = 0;
        if ((x == 1) || (x == LENGTH-2)) if ((j == 1) || (j == 3)) sideflag = 1;
		xp = x+2*dx[j];
		yp = y+2*dy[j];
		
		if (xp > 0 && xp < LENGTH-1 && yp > 0 && yp < WIDTH-1 && lat[xp][yp] == t)
			if (NewRandomInteger < prob || sideflag == 1)  
			{	
				lat[x+dx[j]][y+dy[j]] = s; //activate FK bond for walk to possibly step on
				PutOnStack0(xp,yp)
			}
	}
}

void walk(void)
{
	
	long	x, y, dir, cntr;
	long	end1x, end1y, end2x, end2y;
	long	dx[4] = {-1 ,1, 1, -1}, dy[4] = {1, 1, -1, -1};
	
	dir = 41;
	x = 1;
	y = 0;
	
	walk1[x][y] = runs;
	
	do{	//first walk

		x = x+dx[dir%4];
		y = y+dy[dir%4];

		if (lat[x][y] >= Q && walk1[x][y] == runs) ++dir;  //occupied, already visited by first walk
		else if (lat[x][y] >= Q && walk1[x][y] < runs) { //occupied, first visit by first walk
			PutOnStack1(x,y)
			walk1[x][y] = runs; //visited by first walk
			++dir;
		}
		else if (lat[x][y] < Q && walk1[x][y] == runs) --dir; //vacant, already visited by first walk
		else if (lat[x][y] < Q && walk1[x][y] < runs) { //vacant, first visit by first walk
			PutOnStack1(x,y)
			walk1[x][y] = runs; //visited by first walk
			--dir;
		}
		else {
			printf("error1\n"); 
			exit(1);
		}
		
	}while (((x != 1) || (y != WIDTH-1)) && ((x != LENGTH-2) || (y != 0))); //while we haven't reached the top left or bottom right corner
	end1x = x;
	end1y = y;
	
	dir = 43;
	x = LENGTH-2;
	y = WIDTH-1;
	
	walk2[x][y] = runs;
	
	cntr = 0;
	
	do{	//second walk

		x = x+dx[dir%4];
		y = y+dy[dir%4];

		if (lat[x][y] >= Q && walk2[x][y] == runs) ++dir; //occupied, already visited by second walk (each site can only be "visited" twice)
		else if (lat[x][y] >= Q && walk1[x][y] == runs) { //occupied, already visited by first walk
			PutOnStack2(x,y)
			++pinchpoint1234[x][y];
			++dir;
		}
		else if (lat[x][y] >= Q && walk2[x][y] < runs && walk1[x][y] < runs) { //occupied, first visit by second walk, never visited by first walk
			PutOnStack2(x,y)
			walk2[x][y] = runs; //visited by second walk
			++dir;
		}
		else if (lat[x][y] < Q && walk2[x][y] == runs) --dir; //vacant, already visited by second walk
		else if (lat[x][y] < Q && walk1[x][y] == runs) { //vacant, already visited by first walk
			PutOnStack2(x,y)
			++pinchpoint1234[x][y];
			--dir;
		}
		else if (lat[x][y] < Q && walk2[x][y] < runs && walk1[x][y] < runs) { //vacant, first visit by second walk, never visited by first walk
			PutOnStack2(x,y)
			walk2[x][y] = runs; //visited by second walk
			--dir;
		}
		else {
			printf("error2\n"); 
			exit(1);
		}
		cntr++;
	}while (((x != 1) || (y != WIDTH-1)) && ((x != LENGTH-2) || (y != 0)));
	end2x = x;
	end2y = y;
	
	if (end1x == LENGTH-2 && end1y == 0 && end2x == 1 && end2y == WIDTH-1) connectivity = 1; //horzontal crossing
	else if (end1x == 1 && end1y == WIDTH-1 && end2x == LENGTH-2 && end2y == 0) connectivity = 2; //vertical crossing
	else {
		printf("error3\n");
		exit(1);
	}
}


void randinit(long seed)
{   
    double a, ee = -1 + 1/2147483648.0;
    long i;
    extern long nd, ra[M+1];
    
    a = seed/2147483648.0;
    for (nd = 0; nd <= M; nd++)
    {
        a *= 16807;
        a += ee * (long)(a);
        if (a >= 1) a += ee;
        ra[nd] = (long) (2147483648.0 * a);
    }
    for(i = 1; i<10000001; i++)
        NewRandomInteger;
}