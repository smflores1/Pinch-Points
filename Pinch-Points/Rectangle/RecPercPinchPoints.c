#include "stdlib.h"
#include "stdio.h"
#include "math.h"

#define SEED 811
#define LENGTH 2001 // must be odd
#define WIDTH 1001 //must be odd
#define RUNSMAX 5000005
#define PRINTFREQ 10000
#define PROB 0.5
#define OUTFILE "RecPercPinchPoints.ps"
#define OUTFILE12 "RecPercPinchPoints12.txt"
#define OUTFILE14 "RecPercPinchPoints14.txt"
#define OUTFILE34 "RecPercPinchPoints34.txt"
#define OUTFILE23 "RecPercPinchPoints23.txt"
#define OUTFILE1234 "RecPercPinchPoints1234.txt"
#define RUNFILE "RunFile.txt"
#define S 1048575//32767
#define M  16383					

#define GetFromStack1(X,Y) {X = collide1x[gptr1 & S]; Y = collide1y[gptr1 & S]; ++gptr1;}
#define PutOnStack1(X,Y)	{collide1x[pptr1 & S]=X; collide1y[pptr1 & S]=Y; ++pptr1; if  (pptr1 == gptr1) {printf("error42\n");break;}}

#define GetFromStack2(X,Y) {X = collide2x[gptr2 & S]; Y = collide2y[gptr2 & S]; ++gptr2;}
#define PutOnStack2(X,Y)	{collide2x[pptr2 & S]=X; collide2y[pptr2 & S]=Y; ++pptr2; if  (pptr2 == gptr2) {printf("error43\n");break;}}

#define NewRandomInteger (++nd,ra[nd&M] = ra[(nd-471)&M]^ra[(nd-1586)&M]^ra[(nd-6988)&M]^ra[(nd-9689)&M])	

void randinit(long seed);

long	F = 3;
long 	ra[M+1], nd;
long	lat1[LENGTH][WIDTH];
long	lat2[LENGTH][WIDTH];
long	pinchpoint12[LENGTH][WIDTH];
long	pinchpoint14[LENGTH][WIDTH];
long	pinchpoint34[LENGTH][WIDTH];
long	pinchpoint23[LENGTH][WIDTH];
long	pinchpoint1234[LENGTH][WIDTH];
long 	collide1x[S+1], collide1y[S+1];
long 	collide2x[S+1], collide2y[S+1];

int main(void)
{
	int		connectivity = 0;
	long	centerx, centery;
	long	dx[4] = {-1 ,1, 1, -1}, dy[4] = {1, 1, -1, -1};
	long 	x, y, dir, prob, runs; 
	long	gptr1, pptr1, gptr2, pptr2;
	long	end1x, end1y, end2x, end2y;
	FILE	*fpsketch, *fp12, *fp14, *fp23, *fp34, *fp1234, *fprunfile;
	
	prob = (long) (2147483648.0*PROB);
	randinit(SEED);
	
	centerx = (LENGTH-1)*1.0/2;
	centery = (WIDTH-1)*1.0/2+1;
	
	for (x = 1; x < LENGTH-1; ++x)
		for (y = 1; y < WIDTH-1; ++y){
			lat1[x][y] = 0;
			lat2[x][y] = 0;
			pinchpoint12[x][y] = pinchpoint14[x][y] = pinchpoint34[x][y] = pinchpoint23[x][y] = pinchpoint1234[x][y] = 0;
		}
	
	for(runs = 1; runs <= RUNSMAX; ++runs){
		
		gptr1 = pptr1 = gptr2 = pptr2 = 0;
		
		for (y = 1; y < WIDTH-1; ++y){
			lat1[1][y] = 2*runs-1; //always occupied
			lat2[1][y] = 2*runs-1; //always occupied
			lat1[LENGTH-2][y] = 2*runs-1; //always occupied
			lat2[LENGTH-2][y] = 2*runs-1; //always occupied
		}	
		for (x = 1; x < LENGTH-1; ++x){
			lat1[x][0] = 2*runs; //always vaccant (so first hull walk reflects off "dual wired" bottom row)
			lat2[x][0] = 2*runs; //always vaccant (so second hull walk reflects off "dual wired" bottom row)
			lat1[x][WIDTH-1] = 2*runs; //always vaccant (so first hull walk reflects off "dual wired" top row)
			lat2[x][WIDTH-1] = 2*runs; //always vaccant (so second hull walk reflects off "dual wired" top row)
		}
		
		dir = 41;
		x = 1;
		y = 0;
		
		do{
			
			x = x+dx[dir%4];
			y = y+dy[dir%4];
			
			if (lat1[x][y] > 2*runs) printf("error1\n");
			if (lat1[x][y] == 2*runs) --dir; //already set vacant by first walk
			else if (lat1[x][y] == 2*runs-1) ++dir; //already set occupied by first walk
			else if (NewRandomInteger < prob){
				PutOnStack1(x,y)
				lat1[x][y] = 2*runs-1; //set occupied by first walk
				++dir;
			}
			else {
				PutOnStack1(x,y)
				lat1[x][y] = 2*runs; //set vacant by first walk
				--dir;
			}
			
		}while (((x != 1) || (y != WIDTH-1)) && ((x != LENGTH-2) || (y != 0))); //while we haven't reached the top left or bottom right corner
		end1x = x;
		end1y = y;
		
		dir = 43;
		x = LENGTH-2;
		y = WIDTH-1;
		
		do{
			
			x = x+dx[dir%4];
			y = y+dy[dir%4];
			
			if (lat2[x][y] > 2*runs) printf("error2\n");
			if (lat2[x][y] == 2*runs) --dir; //alread set vacant by second walk
			else if (lat2[x][y] == 2*runs-1) ++dir; //already set occupied by second walk
			else if (lat1[x][y] == 2*runs) { //already set vacant by first walk
				PutOnStack2(x,y)
				lat2[x][y] == 2*runs;
				++pinchpoint1234[x][y];
				--dir;
			}
			else if (lat1[x][y] == 2*runs-1) { //already set occupied by first walk
				PutOnStack2(x,y)
				lat2[x][y] == 2*runs-1;
				++pinchpoint1234[x][y];
				++dir;
			}
			else if (NewRandomInteger < prob){
				PutOnStack2(x,y)
				lat2[x][y] = 2*runs-1; //set occupied by second walk
				++dir;
			}
			else {
				PutOnStack2(x,y)
				lat2[x][y] = 2*runs; //set vacant by second walk
				--dir;
			}
			
		}while (((x != 1) || (y != WIDTH-1)) && ((x != LENGTH-2) || (y != 0)));
		end2x = x;
		end2y = y;
		
		if (end1x == LENGTH-2 && end1y == 0 && end2x == 1 && end2y == WIDTH-1) connectivity = 1; //horzontal crossing
		else if (end1x == 1 && end1y == WIDTH-1 && end2x == LENGTH-2 && end2y == 0) connectivity = 2; //vertical crossing
		else printf("error3\n");
		
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
		
	}/*for runs*/
	
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
	nd = M;
	for(i = 0; i<100001; i++)
		NewRandomInteger;
}


