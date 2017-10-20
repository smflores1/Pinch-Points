#include <fstream>
#include <iostream>
#include <cstdlib>
#include <iomanip>
#include <cmath>

#define length 2002
#define height 2002

#define Tx(X,Y) (X + Y*1.0/2) // transform hexagon shape to a regular hexagon, x coordinate
#define Ty(X,Y) (Y*sqrt(3)*1.0/2) // transform hexagon shape to a regular hexagon, y coordinate

#define Tfullx(X,Y) (Tx(X,Y)/1000 - Tx(length/2,height/2)/1000) //center hexagon, x coordinate
#define Tfully(X,Y) (Ty(X,Y)/1000 - Ty(length/2,height/2)/1000) //center hexagon, y coordinate

int main()
{
	using namespace std;
	
	double first1, second1, third1;
	double centervalue;
	
	ifstream in_stream1;
	ofstream out_stream;
	
	in_stream1.open("HexPercPinchPoints12_34_56.txt"); //change name to your input file
	out_stream.open("HexPercPinchPoints12_34_56Scaled.txt"); //change name to your output file
	
	//find hexagon center and save the value stored at that point
	while (! in_stream1.eof()){
		
		in_stream1 >> first1 >> second1 >> third1;
		
		if (first1 == length/2 && second1 == height/2) 
		{
			centervalue = third1;
			break;
		}
	}
	
	printf("here\n");
	
	in_stream1.close();
	out_stream.close();
	
	in_stream1.open("HexPercPinchPoints12_34_56.txt"); //change name to your input file
	out_stream.open("HexPercPinchPoints12_34_56Scaled.txt"); //change name to your ouput file
	
	while (! in_stream1.eof()){
		
		in_stream1 >> first1 >> second1 >> third1;
		 
		if (third1 != 2147483646 && (int)(first1)%40 == 0  && (int)(second1)%40 == 0) { // trim lattice sites outside hexagon but inside rectangle bounding hexagon
			if (first1 >= 1 && first1 < length){
				if (second1 >= fmax((height-1)/2+1-first1,1) && second1 < fmin(height-1,(height-1)/2+height-1-first1)) {
					out_stream << Tfullx(first1,second1) << "          " << Tfully(first1,second1) << "          " << third1*1.0/centervalue << endl; //transform hexagon so it is regular, normalize lattice variables by center value
				}
			}
		}
	}
	
	in_stream1.close();
	out_stream.close();
	
}
