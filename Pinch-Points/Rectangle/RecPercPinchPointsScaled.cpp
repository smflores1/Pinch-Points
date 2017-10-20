#include <fstream>
#include <iostream>
#include <cstdlib>
#include <iomanip>

int main()
{
	using namespace std;
	
	int cntr = 0;
	
	double first1, second1, third1;
    double centervalue;
	
	ifstream in_stream1;
	
	ofstream out_stream;
	
	in_stream1.open("RecPercPinchPoints12Avg.txt"); //change name to your input file

	out_stream.open("RecPercPinchPoints12Scaled.txt"); //change name to your ouput file
    
    in_stream1 >> first1 >> second1 >> third1;
    
    centervalue = third1;
	
	while (! in_stream1.eof()){
		
		in_stream1 >> first1 >> second1 >> third1;
        
        if (((int) (first1))%2 == 1)
        {

            ++cntr;
    
            if (cntr%20 == 0) out_stream << first1*1.0/1000 << "          " << second1*1.0/1000 << "          " << third1*1.0/centervalue << endl;
            
        }
		
	}
	
	in_stream1.close();

	out_stream.close();
	
}
