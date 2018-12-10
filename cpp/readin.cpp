#include <iostream>
#include <fstream>
using namespace std;
int main(int argc,char* argv[]) {
	string f2o = argv[1];
	string buffer;
	cout<<f2o<<endl;
	ifstream infile; 
	infile.open(f2o); 
	infile>>buffer;
	cout<<buffer<<endl;
	return 0;
}
