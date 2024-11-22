#include <stdlib.h>
#include <vector>
#include <math.h>
#include <iostream>
#include <fstream>
#include <random>
#include <chrono>
#include <string>
#include <iomanip>
#include <sstream>
#include <curand.h>

using namespace std;

typedef double mreal;
typedef vector<mreal> vetor;
#define Ith(v,i) v.at(i)

// blocksxthreads/block = tUnits 

//CaRU Model variabels indice definitions-
#define _numM_          11
#define	_CaRU_numODE_	17
#define _d_             0
#define _f_             1
#define _f2_            2
#define _fCass_         3
#define _R_             4 
#define _O_             5 
#define _RI_            6 
#define _I_             7 
#define _CaSR_          8 
#define _CaSS_          9 
#define _Cai_           10

#define _C_             11
#define _Cl_            12
#define _If_            13
#define _Ifl_           14
#define _If2_           15
#define _If2l_          16
//--------------------------------------------------

class CaRU {
    public:
        int id;
        mreal ICal;
        mreal INaca;
        mreal IbCa;
        mreal IpCa;
        mreal V;
        mreal Nai;
        mreal Cai;
        mreal CaSR;
        __device__ void solveStep(int dt);
    private:
        __device__ void calcDerivatives(int dt);
};

__global__ void solveCaRUStep(CaRU *ca_units, int tUnits, int dt);
