#include <fstream>
#include <iostream>
#include <cstdlib>
#include <iomanip>

int main()
{
	using namespace std;
	
	int cntr = 0;
	
	double first, second, third;
    double centervalue;
	
	ifstream in_stream;
	ofstream out_stream;
	
	in_stream.open("RectFkQ2R2PinchPoints12.txt");
	out_stream.open("RectFkQ2R2PinchPoints12Scaled.txt");
    
    in_stream >> first >> second >> third;
    
    centervalue = third;
	
	while (! in_stream.eof()){
		
		in_stream >> first >> second >> third;
        
        if (((int) (first))%2 == 1)
        {
		
            ++cntr;
		
            if (cntr%20 == 0) out_stream << first*1.0/1000 << "    " << second*1.0/1000 << "    " << third*1.0/centervalue << "    " << endl;
            
        }
		
	}
	
	in_stream.close();

	out_stream.close();
	
}
