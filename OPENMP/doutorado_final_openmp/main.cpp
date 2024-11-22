#include <iostream>
#include <string>
#include "Myocyte.h"

#define	initTime	0.0e0	// ms
#define	endTime		1.0e1	// ms
#define	stepTime	1.0e-3	// ms
#define	printRate	1.0e-4	// 1
#define	saveRate	1.0e-3 	// 1

#define xUnits		16
#define yUnits		16
#define n_LCC		-1
#define n_RyR		n_LCC*5
#define s_LCC		true
#define s_RyR		true
#define model		Myocyte::TM
#define DCai		1.*0.4		// 1/ms
#define DCaSR		1.*0.12	// 1/ms

using namespace::std;

/*
How to run: mpirun -np n_core ./exec n_lcc NxN
*/

int main(int argc, char* argv[]){
	// srandom(time(NULL));
	// srandom(seed);x

	auto begin = chrono::high_resolution_clock::now();
	ios_base::sync_with_stdio(false);

	string arg = argv[1];
	string arg1 = argv[2];
	// string outputPath = "../exit/final/dif0/realizations/"+arg+"_lcc_"+arg1+"/";
	string outputPath = "../exit/"+arg+"_lcc_"+arg1+"/";
	int newNlcc = atoi(argv[1]);
	int newNryr = newNlcc*5;

	Myocyte::setOutputPath(outputPath);

	Myocyte* m = new Myocyte(xUnits,yUnits,newNlcc,newNryr,s_LCC,s_RyR,model,DCai,DCaSR);	

	m->solve(stepTime,initTime,endTime,printRate,saveRate);

	delete m;

	auto end = chrono::high_resolution_clock::now();

	auto duration = chrono::duration_cast<chrono::nanoseconds>(end - begin);

	ofstream file;
    file.open("6_thread(s).txt", ios::app);
    file.precision(8);
	file << "main: " << fixed << (duration.count() * 1e-9)<< "s\n";
	file.close();

    return 0;
}
