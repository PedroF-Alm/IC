#include "../header/CaRU.h"

__device__ void CaRU::solveStep(int dt) {
    this->id *= 2;
    calcDerivatives(dt);
}

__device__ void CaRU::calcDerivatives(int dt) {
    // solveAlgEquations(caru, dt);
    this->id *= 2;
}


__global__ void solveCaRUStep(CaRU *ca_units, int tUnits, int dt) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x + blockIdx.y * blockDim.y;
    
    if (idx >= tUnits)
        return;

    ca_units[idx].id = idx;
    ca_units[idx].solveStep(dt);
    // for(int c=0; c<_CaRU_numODE_; c++) Ith(y,c) = Ith(y,c) + dt*Ith(dy,c);
    // if(d->save){
    
    //     // d->save_y.push_back(d->y);

    // //     mreal curs[8] = {I_up, I_leak, I_rel, I_xfer, I_CaL, I_bCa, I_pCa, I_NaCa};    
    // //     save_c.push_back(curs);

    // //     vetor extras;
    // //     if(!LCC_STO) extras.push_back(yoLCC);
    // //     else extras.push_back(((mreal)oLCCSt)/((mreal)n_LCC));
    // //     if(!RyR_STO) extras.push_back(Ith(y,_O_));
    // //     else extras.push_back(((mreal)oRyRSt)/((mreal)n_RyR));
    // //     extras.push_back(Cai_free);
    // //     extras.push_back(Cass_free);
    // //     extras.push_back(CaSR_free);
    // //     save_e.push_back(extras);
    // }
}