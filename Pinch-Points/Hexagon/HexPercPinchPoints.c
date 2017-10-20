#include "stdlib.h"
#include "stdio.h"
#include "math.h"

void randinit(long seed);

#define SEED 12
#define L1 1000
#define L2 1000
#define L3 1000 
#define L4 1000 
#define RUNSMAX 50000005
#define PRINTFREQ 1000000
#define PROB 0.5
#define OUTFILE "HexPinchPoints.ps"
#define OUTFILE12_34_56 "HexPinchPoints12_34_56.txt"
#define OUTFILE12_36_45 "HexPinchPoints12_36_45.txt"
#define OUTFILE14_23_56 "HexPinchPoints14_23_56.txt"
#define OUTFILE16_23_45 "HexPinchPoints16_23_45.txt"
#define OUTFILE16_25_34 "HexPinchPoints16_25_34.txt"
#define OUTFILE23_14_56 "HexPinchPoints23_14_56.txt"
#define OUTFILE23_16_45 "HexPinchPoints23_16_45.txt"
#define OUTFILE25_16_34 "HexPinchPoints25_16_34.txt"
#define OUTFILE34_12_56 "HexPinchPoints34_12_56.txt"
#define OUTFILE34_16_25 "HexPinchPoints34_16_25.txt"
#define OUTFILE36_12_45 "HexPinchPoints36_12_45.txt"
#define OUTFILE45_12_36 "HexPinchPoints45_12_36.txt"
#define OUTFILE45_16_23 "HexPinchPoints45_16_23.txt"
#define OUTFILE56_12_34 "HexPinchPoints56_12_34.txt"
#define OUTFILE56_14_23 "HexPinchPoints56_14_23.txt"
#define OUTFILE1234 "HexPinchPoints1234.txt"
#define OUTFILE1236 "HexPinchPoints1236.txt"
#define OUTFILE1256 "HexPinchPoints1256.txt"
#define OUTFILE1456 "HexPinchPoints1456.txt"
#define OUTFILE2345 "HexPinchPoints2345.txt"
#define OUTFILE3456 "HexPinchPoints3456.txt"
#define OUTFILE123456 "HexPinchPoints123456.txt"
#define RUNFILE "RunFile.txt"
#define DIRMAX 6
#define S 1048575//32767
#define M  16383					
#define BLOCK 2147483647  //"blocked," we never determine if the site is occupied or vaccant
#define OUTSIDE 2147483646  //"blocked," we never determine if the site is occupied or vaccant

#define GetFromStack1(X,Y) {X = collide1x[gptr1 & S]; Y = collide1y[gptr1 & S]; ++gptr1;}
#define PutOnStack1(X,Y)	{collide1x[pptr1 & S]=X; collide1y[pptr1 & S]=Y; ++pptr1; if  (pptr1 == gptr1) {printf("error42\n");break;}}

#define GetFromStack2(X,Y) {X = collide2x[gptr2 & S]; Y = collide2y[gptr2 & S]; ++gptr2;}
#define PutOnStack2(X,Y)	{collide2x[pptr2 & S]=X; collide2y[pptr2 & S]=Y; ++pptr2; if  (pptr2 == gptr2) {printf("error43\n");break;}}

#define GetFromStack3(X,Y) {X = collide3x[gptr3 & S]; Y = collide3y[gptr3 & S]; ++gptr3;}
#define PutOnStack3(X,Y)	{collide3x[pptr3 & S]=X; collide3y[pptr3 & S]=Y; ++pptr3; if  (pptr3 == gptr3) {printf("error44\n");break;}}

#define GetFromStack12(X,Y) {X = collide12x[gptr12 & S]; Y = collide12y[gptr12 & S]; ++gptr12;}
#define PutOnStack12(X,Y)	{collide12x[pptr12 & S]=X; collide12y[pptr12 & S]=Y; ++pptr12; if  (pptr12 == gptr12) {printf("error45\n");break;}}

#define GetFromStack13(X,Y) {X = collide13x[gptr13 & S]; Y = collide13y[gptr13 & S]; ++gptr13;}
#define PutOnStack13(X,Y)	{collide13x[pptr13 & S]=X; collide13y[pptr13 & S]=Y; ++pptr13; if  (pptr13 == gptr13) {printf("error46\n");break;}}

#define GetFromStack23(X,Y) {X = collide23x[gptr23 & S]; Y = collide23y[gptr23 & S]; ++gptr23;}
#define PutOnStack23(X,Y)	{collide23x[pptr23 & S]=X; collide23y[pptr23 & S]=Y; ++pptr23; if  (pptr23 == gptr23) {printf("error47\n");break;}}

#define Tx(X,Y) (X + Y*1.0/2)
#define Ty(X,Y) (Y*sqrt(3)*1.0/2)

#define NewRandomInteger (++nd,ra[nd&M] = ra[(nd-471)&M]^ra[(nd-1586)&M]^ra[(nd-6988)&M]^ra[(nd-9689)&M])								 

long	F = 4;
long 	ra[M+1], nd;
long	L5 = L1+L2-L4;
long	L6 = L3+L4-L1;
long	vertex1x = L1+1;
long	vertex1y = 1;
long	vertex2x = L1+L2+1;
long	vertex2y = 1;
long	vertex3x = L1+L2+1;
long	vertex3y = L3+1;
long	vertex4x = L1+L2-L4+1;
long	vertex4y = L3+L4+1;
long	vertex5x = 1;
long	vertex5y = L3+L4+1;
long	vertex6x = 1;
long	vertex6y = L1+1;
long	length = L1+L2+2;
long	height = L3+L4+2;
long	lat1[L1+L2+3][L3+L4+3];
long	lat2[L1+L2+3][L3+L4+3];
long	lat3[L1+L2+3][L3+L4+3];
long 	collide1x[S+1], collide1y[S+1];
long 	collide2x[S+1], collide2y[S+1];
long 	collide3x[S+1], collide3y[S+1];
long 	collide12x[S+1], collide12y[S+1];
long 	collide13x[S+1], collide13y[S+1];
long 	collide23x[S+1], collide23y[S+1];
long	pinchpoint12_34_56[L1+L2+3][L3+L4+3];
long	pinchpoint12_36_45[L1+L2+3][L3+L4+3];
long	pinchpoint14_23_56[L1+L2+3][L3+L4+3];
long	pinchpoint16_23_45[L1+L2+3][L3+L4+3];
long	pinchpoint16_25_34[L1+L2+3][L3+L4+3];
long	pinchpoint23_14_56[L1+L2+3][L3+L4+3];
long	pinchpoint23_16_45[L1+L2+3][L3+L4+3];
long	pinchpoint25_16_34[L1+L2+3][L3+L4+3];
long	pinchpoint34_12_56[L1+L2+3][L3+L4+3];
long	pinchpoint34_16_25[L1+L2+3][L3+L4+3];
long	pinchpoint36_12_45[L1+L2+3][L3+L4+3];
long	pinchpoint45_12_36[L1+L2+3][L3+L4+3];
long	pinchpoint45_16_23[L1+L2+3][L3+L4+3];
long	pinchpoint56_12_34[L1+L2+3][L3+L4+3];
long	pinchpoint56_14_23[L1+L2+3][L3+L4+3];
long	pinchpoint1234[L1+L2+3][L3+L4+3];
long	pinchpoint1236[L1+L2+3][L3+L4+3];
long	pinchpoint1256[L1+L2+3][L3+L4+3];
long	pinchpoint1456[L1+L2+3][L3+L4+3];
long	pinchpoint2345[L1+L2+3][L3+L4+3];
long	pinchpoint3456[L1+L2+3][L3+L4+3];
long	pinchpoint123456[L1+L2+3][L3+L4+3];


int main(void)

{
	int connectivity = 0;
	
	long	dx[6] = {-1, -1, 0, 1, 1, 0}, dy[6] = {0, 1, 1, 0, -1, -1};
	long 	x, y, i, xp, yp, dir, prob;
	long	gptr1, pptr1, gptr2, pptr2, gptr3, pptr3, gptr12, pptr12, gptr13, pptr13, gptr23, pptr23; 
	long	end1x, end1y, end2x, end2y, end3x, end3y;
	long	runs;
	double	 x1, x2, y1, y2;
	
	FILE    *fp12_34_56, *fp12_36_45, *fp14_23_56, *fp16_23_45, *fp16_25_34, *fp23_14_56, *fp23_16_45, *fp25_16_34, *fp34_12_56, *fp34_16_25;
	FILE	*fp36_12_45, *fp45_12_36, *fp45_16_23, *fp56_12_34, *fp56_14_23, *fp1234, *fp1236, *fp1256, *fp1456, *fp2345, *fp3456, *fp123456;
	FILE	*fprunfile;
	
	prob = (long) (2147483648.0*PROB);
	randinit(SEED);
	
	if (L4 > L1+L2) {
		printf("error40\n");
		exit(1);
	}
	
	if (L1 > L3+L4) {
		printf("error41\n");
		exit(1);
	}
	
	for (x = 0; x <= length; ++x)
		for (y = 0; y <= height; ++y)		
		{
			
			lat1[x][y] = lat2[x][y] = lat3[x][y] = 0;
			pinchpoint12_34_56[x][y] = pinchpoint12_36_45[x][y] = pinchpoint14_23_56[x][y] = 0;
			pinchpoint16_23_45[x][y] = pinchpoint16_25_34[x][y] = pinchpoint23_14_56[x][y] = 0;
			pinchpoint23_16_45[x][y] = pinchpoint25_16_34[x][y] = pinchpoint34_12_56[x][y] = 0;
			pinchpoint34_16_25[x][y] = pinchpoint36_12_45[x][y] = pinchpoint45_12_36[x][y] = 0;
			pinchpoint45_16_23[x][y] = pinchpoint56_12_34[x][y] = pinchpoint56_14_23[x][y] = 0;
			pinchpoint1234[x][y] = pinchpoint1236[x][y] = pinchpoint1256[x][y] = 0;
			pinchpoint1456[x][y] = pinchpoint2345[x][y] = pinchpoint3456[x][y] = 0;
			pinchpoint123456[x][y] = 0;
			
			if (x+y == vertex1x){
				lat1[x][y] = lat2[x][y] = lat3[x][y] = OUTSIDE; //always blocked
			}
			if ((x+y == vertex1x+1) && (x > 1) && (y > 0)){
				lat1[x][y] = lat2[x][y] = lat3[x][y] = BLOCK; //always vacant
			}
			
			if ((vertex1x <= x) && (y == 0)){
				lat1[x][y] = lat2[x][y] = lat3[x][y] = OUTSIDE; //always blocked
			}
			
			if ((y <= vertex3y) && (x == length)){
				lat1[x][y] = lat2[x][y] = lat3[x][y] = OUTSIDE; //always blocked
			}
			if ((y <= vertex3y) && (x == length-1) && (y > 1)){
				lat1[x][y] = lat2[x][y] = lat3[x][y] = BLOCK; //always vacant
			}
			
			if (x+y == length+L3+1){
				lat1[x][y] = lat2[x][y] = lat3[x][y] = OUTSIDE; //always blocked
			}
			
			if ((x <= vertex4x) && (y == height)){ 
				lat1[x][y] = lat2[x][y] = lat3[x][y] = OUTSIDE; //alwyas blocked
			}
			if ((x > 0) && (x < vertex4x) && (y == height-1)){
				lat1[x][y] = lat2[x][y] = lat3[x][y] = BLOCK; //always vacant
			}
			
			if ((vertex6y <= y) && (y <= height) && (x == 0)){
				lat1[x][y] = lat2[x][y] = lat3[x][y] = OUTSIDE; //always blocked
			}
			
		}
	
	for (runs = 1; runs <= RUNSMAX; ++runs)
	{
		
		gptr1 = pptr1 = gptr2 = pptr2 = gptr3 = pptr3 = gptr12 = pptr12 = gptr13 = pptr13 = gptr23 = pptr23 = 0;
		
		
		for (x = vertex1x+1; x < length; ++x){
			lat1[x][1] = 2*runs-1; //occupied
		}
		
		for (y = vertex6y; y < vertex5y; ++y){
			lat2[1][y] = 2*runs-1; //occupied
		}
		
		for (x = vertex4x; x < length-1; ++x){
			y = length+L3-x;
			lat3[x][y] = 2*runs-1; //occupied
		}
		
		dir = 60;
		x = vertex1x+1;
		y = 1;
		
		do //walk 1
		{
			xp = x+dx[dir%6];
			yp = y+dy[dir%6];
			if ((lat1[xp][yp] >= 2*runs) || (lat2[xp][yp] >= 2*runs) || (lat3[xp][yp] >= 2*runs)) { //either blocked or already set vacant
				if (lat1[xp][yp] != BLOCK && lat1[xp][yp] != OUTSIDE) lat1[xp][yp] = 2*runs;
				++dir;
			}
			else if ((lat1[xp][yp] == 2*runs-1) || (lat2[xp][yp] == 2*runs-1) || (lat3[xp][yp] == 2*runs-1) || NewRandomInteger < prob) { //occupied
				PutOnStack1(xp,yp) //a one pinch point
				lat1[xp][yp] = 2*runs-1; //occupied
				x = xp; 
				y = yp;
				--dir;
			}
			else { //set vacant by first walk
				lat1[xp][yp] = 2*runs;
				++dir;
			}
		}
		while (((x != vertex2x) || (y != vertex2y)) && ((x != vertex4x) || (y != vertex4y)) && ((x != vertex6x) || (y != vertex6y)));
		end1x = x;
		end1y = y;
		
		dir = 62;
		x = vertex5x;
		y = vertex5y-1;
		
		do //walk 2
		{
			xp = x+dx[dir%6];
			yp = y+dy[dir%6];
			if ((lat1[xp][yp] >= 2*runs) || (lat2[xp][yp] >= 2*runs) || (lat3[xp][yp] >= 2*runs)) { //either blocked or already set vacant
				if (lat1[xp][yp] == 2*runs) { //a two pinch point between first and second walks (deflection)
					PutOnStack12(xp,yp)
					PutOnStack2(xp,yp) 
				}
				if (lat2[xp][yp] != BLOCK && lat2[xp][yp] != OUTSIDE) lat2[xp][yp] = 2*runs;
				++dir;
			}
			else if ((lat1[xp][yp] == 2*runs-1) || (lat2[xp][yp] == 2*runs-1) || (lat3[xp][yp] == 2*runs-1) || NewRandomInteger < prob) { //occupied
				if (lat1[xp][yp] == 2*runs-1){ //a two pinch point between first and second walks (merging)
					PutOnStack12(xp,yp)
					PutOnStack2(xp,yp) 
				}
				else { //a one pinch point
					PutOnStack2(xp,yp)
				}
				lat2[xp][yp] = 2*runs-1; //occupied
				x = xp;
				y = yp;
				--dir;
			}
			else { //set vacant by second walk
				lat2[xp][yp] = 2*runs;
				++dir;
			}
		}
		while (((x != vertex2x) || (y != vertex2y)) && ((x != vertex4x) || (y != vertex4y)) && ((x != vertex6x) || (y != vertex6y)));
		end2x = x;
		end2y = y;
		
		dir = 64;
		x = vertex2x-1;
		y = vertex3y+1;
		
		do //walk 3
		{
			xp = x+dx[dir%6];
			yp = y+dy[dir%6];
			if ((lat1[xp][yp] >= 2*runs) || (lat2[xp][yp] >= 2*runs) || (lat3[xp][yp] >= 2*runs)) { //either blocked or already set vacant
				if (lat1[xp][yp] == 2*runs && lat2[xp][yp] != 2*runs) { //a two pinch point between first and third walks (deflection)
					PutOnStack13(xp,yp)
					PutOnStack3(xp,yp)
				}
				else if (lat2[xp][yp] == 2*runs && lat1[xp][yp] != 2*runs) { //a two pinch point between second and third walks (deflection)
					PutOnStack23(xp,yp)
					PutOnStack3(xp,yp)
				}
				else if (lat1[xp][yp] == 2*runs && lat2[xp][yp] == 2*runs) { //a three pinch point between first and second and third walks (deflection)
					++pinchpoint123456[xp][yp];
					PutOnStack13(xp,yp)
					PutOnStack23(xp,yp)
					PutOnStack3(xp,yp)
				}
				if (lat3[xp][yp] != BLOCK && lat3[xp][yp] != OUTSIDE) lat3[xp][yp] = 2*runs;
				++dir;
			}
			else if ((lat1[xp][yp] == 2*runs-1) || (lat2[xp][yp] == 2*runs-1) || (lat3[xp][yp] == 2*runs-1) || NewRandomInteger < prob) { //occupied
				
				lat3[xp][yp] = 2*runs-1; //occupied
				if (lat1[xp][yp] == 2*runs-1 && lat2[xp][yp] != 2*runs-1  && yp != 1){ //a two pinch point between first and third walks (merging)
					PutOnStack13(xp,yp)
					PutOnStack3(xp,yp)
				}
				else if (lat2[xp][yp] == 2*runs-1 && lat1[xp][yp] != 2*runs-1  && xp != 1){ //a two pinch point between second and third walks (merging)
					PutOnStack23(xp,yp)
					PutOnStack3(xp,yp)
				}
				else if (lat1[xp][yp] == 2*runs-1 && lat2[xp][yp] == 2*runs-1){ //a three pinch point between first and second and third walks (merging)
					++pinchpoint123456[xp][yp];
					PutOnStack13(xp,yp)
					PutOnStack23(xp,yp)
					PutOnStack3(xp,yp)
				}
				else { //a one pinch point
					PutOnStack3(xp,yp) //a one pinch point
				}
				x = xp;
				y = yp;
				--dir;				
			}
			else { //set vacant by third walk
				lat3[xp][yp] = 2*runs;
				++dir;
			}
		}
		while (((x != vertex2x) || (y != vertex2y)) && ((x != vertex4x) || (y != vertex4y)) && ((x != vertex6x) || (y != vertex6y)));
		end3x = x;
		end3y = y;
		
		if (end1x == vertex2x && end1y == vertex2y && end2x == vertex6x && end2y == vertex6y && end3x == vertex4x && end3y == vertex4y) connectivity = 1;
		else if (end1x == vertex6x && end1y == vertex6y && end2x == vertex4x && end2y == vertex4y && end3x == vertex2x && end3y == vertex2y) connectivity = 2;
		else if (end1x == vertex4x && end1y == vertex4y && end2x == vertex6x && end2y == vertex6y && end3x == vertex2x && end3y == vertex2y) connectivity = 3;
		else if (end1x == vertex6x && end1y == vertex6y && end2x == vertex2x && end2y == vertex2y && end3x == vertex4x && end3y == vertex4y) connectivity = 4;
		else if (end1x == vertex2x && end1y == vertex2y && end2x == vertex4x && end2y == vertex4y && end3x == vertex6x && end3y == vertex6y) connectivity = 5;
		
		while (gptr1 != pptr1) {
			GetFromStack1(x,y)
			if (connectivity == 1) ++pinchpoint12_34_56[x][y];
			if (connectivity == 2) ++pinchpoint16_23_45[x][y];
			if (connectivity == 3) ++pinchpoint14_23_56[x][y];
			if (connectivity == 4) ++pinchpoint16_25_34[x][y];
			if (connectivity == 5) ++pinchpoint12_36_45[x][y];
		}
		
		while (gptr2 != pptr2) {
			GetFromStack2(x,y)
			if (connectivity == 1) ++pinchpoint56_12_34[x][y];
			if (connectivity == 2) ++pinchpoint45_16_23[x][y];
			if (connectivity == 3) ++pinchpoint56_14_23[x][y];
			if (connectivity == 4) ++pinchpoint25_16_34[x][y];
			if (connectivity == 5) ++pinchpoint45_12_36[x][y];
		}
		
		while (gptr3 != pptr3) {
			GetFromStack3(x,y)
			if (connectivity == 1) ++pinchpoint34_12_56[x][y];
			if (connectivity == 2) ++pinchpoint23_16_45[x][y];
			if (connectivity == 3) ++pinchpoint23_14_56[x][y];
			if (connectivity == 4) ++pinchpoint34_16_25[x][y];
			if (connectivity == 5) ++pinchpoint36_12_45[x][y];
		}
		
		while (gptr12 != pptr12) {
			GetFromStack12(x,y)
			if (connectivity == 1) ++pinchpoint1256[x][y];
			if (connectivity == 2) ++pinchpoint1456[x][y];
			if (connectivity == 3) ++pinchpoint1456[x][y];
			if (connectivity == 4) ++pinchpoint1256[x][y];
		}
		
		while (gptr13 != pptr13) {
			GetFromStack13(x,y)
			if (connectivity == 1) ++pinchpoint1234[x][y];
			if (connectivity == 2) ++pinchpoint1236[x][y];
			if (connectivity == 3) ++pinchpoint1234[x][y];
			if (connectivity == 5) ++pinchpoint1236[x][y];
		}
		
		while (gptr23 != pptr23) {
			GetFromStack23(x,y)
			if (connectivity == 1) ++pinchpoint3456[x][y];
			if (connectivity == 2) ++pinchpoint2345[x][y];
			if (connectivity == 4) ++pinchpoint2345[x][y];
			if (connectivity == 5) ++pinchpoint3456[x][y];
		}
				
		if (runs%PRINTFREQ == 0)  {
			
			fprunfile = fopen(RUNFILE, "w");
			fprintf(fprunfile, "%15ld\n", runs);
			fclose(fprunfile);
			
			fp12_34_56 = fopen(OUTFILE12_34_56,"w");
			fp12_36_45 = fopen(OUTFILE12_36_45,"w");
			fp14_23_56 = fopen(OUTFILE14_23_56,"w");
			fp16_23_45 = fopen(OUTFILE16_23_45,"w");
			fp16_25_34 = fopen(OUTFILE16_25_34,"w");
			fp23_14_56 = fopen(OUTFILE23_14_56,"w");
			fp23_16_45 = fopen(OUTFILE23_16_45,"w");
			fp25_16_34 = fopen(OUTFILE25_16_34,"w");
			fp34_12_56 = fopen(OUTFILE34_12_56,"w");
			fp34_16_25 = fopen(OUTFILE34_16_25,"w");
			fp36_12_45 = fopen(OUTFILE36_12_45,"w");
			fp45_12_36 = fopen(OUTFILE45_12_36,"w");
			fp45_16_23 = fopen(OUTFILE45_16_23,"w");
			fp56_12_34 = fopen(OUTFILE56_12_34,"w");
			fp56_14_23 = fopen(OUTFILE56_14_23,"w");
			fp1234 = fopen(OUTFILE1234,"w");
			fp1236 = fopen(OUTFILE1236,"w");
			fp1256 = fopen(OUTFILE1256,"w");
			fp1456 = fopen(OUTFILE1456,"w");
			fp2345 = fopen(OUTFILE2345,"w");
			fp3456 = fopen(OUTFILE3456,"w");
			fp123456 = fopen(OUTFILE123456,"w");

			for (x = 0; x <= length; x++)
				for (y = 0; y <= height; y++) {
					fprintf(fp12_34_56, "%15ld%15ld%15ld\n", x, y, pinchpoint12_34_56[x][y]);
					fprintf(fp12_36_45, "%15ld%15ld%15ld\n", x, y, pinchpoint12_36_45[x][y]);
					fprintf(fp14_23_56, "%15ld%15ld%15ld\n", x, y, pinchpoint14_23_56[x][y]);
					fprintf(fp16_23_45, "%15ld%15ld%15ld\n", x, y, pinchpoint16_23_45[x][y]);
					fprintf(fp16_25_34, "%15ld%15ld%15ld\n", x, y, pinchpoint16_25_34[x][y]);
					fprintf(fp23_14_56, "%15ld%15ld%15ld\n", x, y, pinchpoint23_14_56[x][y]);
					fprintf(fp23_16_45, "%15ld%15ld%15ld\n", x, y, pinchpoint23_16_45[x][y]);
					fprintf(fp25_16_34, "%15ld%15ld%15ld\n", x, y, pinchpoint25_16_34[x][y]);
					fprintf(fp34_12_56, "%15ld%15ld%15ld\n", x, y, pinchpoint34_12_56[x][y]);
					fprintf(fp34_16_25, "%15ld%15ld%15ld\n", x, y, pinchpoint34_16_25[x][y]);
					fprintf(fp36_12_45, "%15ld%15ld%15ld\n", x, y, pinchpoint36_12_45[x][y]);
					fprintf(fp45_12_36, "%15ld%15ld%15ld\n", x, y, pinchpoint45_12_36[x][y]);
					fprintf(fp45_16_23, "%15ld%15ld%15ld\n", x, y, pinchpoint45_16_23[x][y]);
					fprintf(fp56_12_34, "%15ld%15ld%15ld\n", x, y, pinchpoint56_12_34[x][y]);
					fprintf(fp56_14_23, "%15ld%15ld%15ld\n", x, y, pinchpoint56_14_23[x][y]);
					fprintf(fp1234, "%15ld%15ld%15ld\n", x, y, pinchpoint1234[x][y]);
					fprintf(fp1236, "%15ld%15ld%15ld\n", x, y, pinchpoint1236[x][y]);
					fprintf(fp1256, "%15ld%15ld%15ld\n", x, y, pinchpoint1256[x][y]);
					fprintf(fp1456, "%15ld%15ld%15ld\n", x, y, pinchpoint1456[x][y]);
					fprintf(fp2345, "%15ld%15ld%15ld\n", x, y, pinchpoint2345[x][y]);
					fprintf(fp3456, "%15ld%15ld%15ld\n", x, y, pinchpoint3456[x][y]);
					fprintf(fp123456, "%15ld%15ld%15ld\n", x, y, pinchpoint123456[x][y]);
				}
			
			fclose(fp12_34_56);
			fclose(fp12_36_45);
			fclose(fp14_23_56);
			fclose(fp16_23_45);
			fclose(fp16_25_34);
			fclose(fp23_14_56);
			fclose(fp23_16_45);
			fclose(fp25_16_34);
			fclose(fp34_12_56);
			fclose(fp34_16_25);
			fclose(fp36_12_45);
			fclose(fp45_12_36);
			fclose(fp45_16_23);
			fclose(fp56_12_34);
			fclose(fp56_14_23);
			fclose(fp1234);
			fclose(fp1236);
			fclose(fp1256);
			fclose(fp1456);
			fclose(fp2345);
			fclose(fp3456);
			fclose(fp123456);
		}
		
	}// for runs
	
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



