#ifndef CARU_H
#define CARU_H

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

//Type definitions--------
typedef double mreal;
typedef vector<mreal> vetor;
#define Ith(v,i) v.at(i)
//------------------------

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

typedef struct {
    mreal y[_CaRU_numODE_];
    mreal dy[_CaRU_numODE_];
    mreal mult[_numM_];
    mreal K_buf_c;
    mreal Buf_c;
    mreal K_buf_sr;
    mreal Buf_sr;
    mreal K_buf_ss;
    mreal Buf_ss;
    mreal V;
    mreal Nai;
    mreal Cai_free;
    mreal CaSR_free;    
    mreal Cass_free;
    mreal RTONF;
    mreal Ca_o;
    mreal tau_d;
    mreal tau_f;
    mreal tau_f2;
    mreal tau_fCass;
    mreal d_inf;
    mreal f_inf;
    mreal f2_inf;
    mreal fCass_inf;
    mreal yoLCC;
    mreal Vmax_up;
    mreal K_up;
    mreal V_leak;
    mreal V_rel;
    mreal V_xfer;
    mreal g_CaL;
    mreal F;
    mreal R;
    mreal T;
    mreal g_bCa;
    mreal g_pCa;
    mreal K_pCa;
    mreal K_NaCa;
    mreal Km_Nai;
    mreal Na_o;
    mreal Km_Ca;
    mreal K_sat;
    mreal n;
    mreal I_up;
    mreal I_leak;
    mreal I_rel;
    mreal I_xfer;
    mreal I_CaL;
    mreal I_bCa;
    mreal I_pCa;
    mreal I_NaCa;
    mreal V_c;
    mreal V_ss;
    mreal V_sr;
    mreal inversevssF2;
    mreal inverseVcF2;
    mreal Cm;
    mreal r_C_If;
    mreal r_C_If2;
    mreal r_C_Cl;
    mreal r_Cl_C;
    mreal r_Cl_Ifl;
    mreal r_Cl_If2l;
    mreal r_Cl_O;
    mreal r_If_C;
    mreal r_If_Ifl;
    mreal r_If2_C;
    mreal r_Ifl_Cl;
    mreal r_Ifl_If;
    mreal r_If2l_Cl;
    mreal r_O_Cl;
    mreal r_O_Ifl;
    mreal r_O_If2l;
    mreal r_Ifl_O;
    mreal r_If2l_If2;
    mreal r_If2_If2l;
    mreal r_If2l_O;
    mreal r_O_R;
    mreal r_RI_R;
    mreal r_R_O;
    mreal r_R_RI;
    mreal r_I_O;
    mreal r_O_I;
    mreal r_RI_I;
    mreal r_I_RI;
    mreal k_NaCa;
    mreal max_sr;
    mreal min_sr;
    mreal EC;
    mreal k1_prime;
    mreal k2_prime;
    mreal k3;
    mreal k4;
    int TT;
    int TM;
    int oLCCSt;
    int oRyRSt;
    int n_LCC;
    int n_RyR;
    int MOD;
    int rRyRSt;
    int riRyRSt;
    int iRyRSt;
    int cLCCSt;
    int ifLCCSt;
    int if2LCCSt;
    int clLCCSt;
    int iflLCCSt;
    int if2lLCCSt;
    bool LCC_STO;
    bool RyR_STO;
    bool save;
} DataBlock;

typedef struct {
    DataBlock *device_c_units;
    mreal *device_V;
    mreal *device_Nai;  
    bool *device_toSave; 
    mreal V; 
    mreal Nai; 
    bool toSave; 
    int tUnits; 
    int dt; int t; 
    int saveRate;
} SolveStepParameters;

class CaRU;

__device__ void updateRyRs(DataBlock *d, mreal dt);
__device__ void updateLCCs(DataBlock *d, mreal dt);
__device__ int calcBinomial_4(int t, mreal p);
__device__ int calcBinomial(int t, mreal p);
__device__ mreal solveQuadratic(mreal a, mreal b, mreal c);
__device__ void solveAlgEquations(DataBlock *d, mreal dt);
__device__ void calcDerivatives(DataBlock *d, mreal dt);
__device__ void solveStep(DataBlock *d, mreal dt);
__global__ void solveCaUnitsStep(SolveStepParameters ssp);

class CaRU{
public:
    CaRU(int id);
    void saveVariables(vetor save_t, string output_path);
    void saveCurrents(vetor save_t, string output_path);
    void saveExtras(vetor save_t, string output_path);
    vetor getVariables();
    vetor getDerivatives();
    mreal getICaL();
    mreal getIbCa();
    mreal getINaCa();
    mreal getIpCa();
    mreal getCaiFree();
    mreal getCaSRFree();
    mreal getCai();
    mreal getCaSR();
    void setCai(mreal value);
    void setCaSR(mreal value);
    void printId();

    //----Model Definitions--   
    const static int TT = 0;
    const static int TM = 1;
    //-----------------------             

    DataBlock toDataBlock() {
        DataBlock d = {};
        for (int i = 0; i < _CaRU_numODE_; i++)
        {
            d.y[i] = this->y[i];
            d.dy[i] = this->dy[i];
        }        
        for (int i = 0; i < _numM_; i++)
        {
            d.mult[i] = this->mult[i];
        }
        d.K_buf_c = this->K_buf_c;
        d.Buf_c = this->Buf_c;
        d.K_buf_sr = this->K_buf_sr;
        d.Buf_sr = this->Buf_sr;
        d.K_buf_ss = this->K_buf_ss;
        d.Buf_ss = this->Buf_ss;
        d.V = this->V;
        d.Nai = this->Nai;
        d.Cai_free = this->Cai_free;
        d.CaSR_free = this->CaSR_free;    
        d.Cass_free = this->Cass_free;
        d.RTONF = this->RTONF;
        d.Ca_o = this->Ca_o;
        d.tau_d = this->tau_d;
        d.tau_f = this->tau_f;
        d.tau_f2 = this->tau_f2;
        d.tau_fCass = this->tau_fCass;
        d.d_inf = this->d_inf;
        d.f_inf = this->f_inf;
        d.f2_inf = this->f2_inf;
        d.fCass_inf = this->fCass_inf;
        d.yoLCC = this->yoLCC;
        d.Vmax_up = this->yoLCC;
        d.K_up = this->K_up;
        d.V_leak = this->V_leak;
        d.V_rel = this->V_rel;
        d.V_xfer = this->V_xfer;
        d.g_CaL = this->g_CaL;
        d.F = this->F;
        d.R = this->R;
        d.T = this->T;
        d.g_bCa = this->g_bCa;
        d.g_pCa = this->g_pCa;
        d.K_pCa = this->K_pCa;
        d.K_NaCa = this->k_NaCa;
        d.Km_Nai = this->Km_Nai;
        d.Na_o = this->Na_o;
        d.Km_Ca = this->Km_Ca;
        d.K_sat = this->K_sat;
        d.n = this->n;
        d.I_up = this->I_up;
        d.I_leak = this->I_leak;
        d.I_rel = this->I_rel;
        d.I_xfer = this->I_xfer;
        d.I_CaL = this->I_CaL;
        d.I_bCa = this->I_bCa;
        d.I_pCa = this->I_pCa;
        d.I_NaCa = this->I_NaCa;
        d.V_c = this->V_c;
        d.V_ss = this->V_ss;
        d.V_sr = this->V_sr;
        d.inversevssF2 = this->inversevssF2;
        d.inverseVcF2 = this->inverseVcF2;
        d.Cm = this->Cm;
        d.r_C_If = this->r_C_If;
        d.r_C_If2 = this->r_C_If2;
        d.r_C_Cl = this->r_C_Cl;
        d.r_Cl_C = this->r_Cl_C;
        d.r_Cl_Ifl = this->r_Cl_Ifl;
        d.r_Cl_If2l = this->r_Cl_If2l;
        d.r_Cl_O = this->r_Cl_O;
        d.r_If_C = this->r_If_C;
        d.r_If_Ifl = this->r_If_Ifl;
        d.r_If2_C = this->r_If2_C;
        d.r_Ifl_Cl = this->r_Ifl_Cl;
        d.r_Ifl_If = this->r_Ifl_If;
        d.r_If2l_Cl = this->r_If2l_Cl;
        d.r_O_Cl = this->r_O_Cl;
        d.r_O_Ifl = this->r_O_Ifl;
        d.r_O_If2l = this->r_O_If2l;
        d.r_Ifl_O = this->r_Ifl_O;
        d.r_If2l_If2 = this->r_If2l_If2;
        d.r_If2_If2l = this->r_If2_If2l;
        d.r_If2l_O = this->r_If2l_O;
        d.r_O_R = this->r_O_R;
        d.r_RI_R = this->r_RI_R;
        d.r_R_O = this->r_R_O;
        d.r_R_RI = this->r_R_RI;
        d.r_I_O = this->r_I_O;
        d.r_O_I = this->r_O_I;
        d.r_RI_I = this->r_RI_I;
        d.r_I_RI = this->r_I_RI;
        d.k_NaCa = this->k_NaCa;
        d.TT = this->TT;
        d.TM = this->TM;
        d.oLCCSt = this->oLCCSt;
        d.oRyRSt = this->oRyRSt;
        d.n_LCC = this->n_LCC;
        d.n_RyR = this->n_RyR;
        d.MOD = this->MOD;
        d.LCC_STO = this->LCC_STO;
        d.RyR_STO = this->RyR_STO;
        d.max_sr = this->max_sr;
        d.min_sr = this->min_sr;
        d.EC = this->EC;
        d.k1_prime = this->k1_prime;
        d.k2_prime = this->k2_prime;
        d.rRyRSt = this->rRyRSt;
        d.riRyRSt = this->riRyRSt;
        d.iRyRSt = this->iRyRSt;
        d.k3 = this->k3;
        d.k4 = this->k4;
        d.cLCCSt = this->cLCCSt;
        d.ifLCCSt = this->ifLCCSt;
        d.if2LCCSt = this->if2LCCSt;
        d.clLCCSt = this->clLCCSt;
        d.iflLCCSt = this->iflLCCSt;
        d.if2lLCCSt = this->if2lLCCSt;
        d.save = this->save;
        return d;
    }

    void updateCaruFromDataBlock(DataBlock *d){
        for (int i = 0; i < _CaRU_numODE_; i++)
        {
            this->y[i] = d->y[i];
            this->dy[i] = d->dy[i];
        }
        this->V = d->V;
        this->Nai = d->Nai;
        this->Cai_free = d->Cai_free;
        this->CaSR_free = d->CaSR_free;    
        this->Cass_free = d->Cass_free;
        this->tau_d = d->tau_d;
        this->tau_f = d->tau_f;
        this->tau_f2 = d->tau_f2;
        this->tau_fCass = d->tau_fCass;
        this->d_inf = d->d_inf;
        this->f_inf = d->f_inf;
        this->f2_inf = d->f2_inf;
        this->fCass_inf = d->fCass_inf;
        this->yoLCC = d->yoLCC;
        this->I_up = d->I_up;
        this->I_leak = d->I_leak;
        this->I_rel = d->I_rel;
        this->I_xfer = d->I_xfer;
        this->I_CaL = d->I_CaL;
        this->I_bCa = d->I_bCa;
        this->I_pCa = d->I_pCa;
        this->I_NaCa = d->I_NaCa;
        this->r_C_If = d->r_C_If;
        this->r_C_If2 = d->r_C_If2;
        this->r_C_Cl = d->r_C_Cl;
        this->r_Cl_C = d->r_Cl_C;
        this->r_Cl_Ifl = d->r_Cl_Ifl;
        this->r_Cl_If2l = d->r_Cl_If2l;
        this->r_Cl_O = d->r_Cl_O;
        this->r_If_C = d->r_If_C;
        this->r_If_Ifl = d->r_If_Ifl;
        this->r_If2_C = d->r_If2_C;
        this->r_Ifl_Cl = d->r_Ifl_Cl;
        this->r_Ifl_If = d->r_Ifl_If;
        this->r_If2l_Cl = d->r_If2l_Cl;
        this->r_O_Cl = d->r_O_Cl;
        this->r_O_Ifl = d->r_O_Ifl;
        this->r_O_If2l = d->r_O_If2l;
        this->r_Ifl_O = d->r_Ifl_O;
        this->r_If2l_If2 = d->r_If2l_If2;
        this->r_If2_If2l = d->r_If2_If2l;
        this->r_If2l_O = d->r_If2l_O;
        this->r_O_R = d->r_O_R;
        this->r_RI_R = d->r_RI_R;
        this->r_R_O = d->r_R_O;
        this->r_R_RI = d->r_R_RI;
        this->r_I_O = d->r_I_O;
        this->r_O_I = d->r_O_I;
        this->r_RI_I = d->r_RI_I;
        this->r_I_RI = d->r_I_RI;
        this->oLCCSt = d->oLCCSt;
        this->oRyRSt = d->oRyRSt;
        this->n_LCC = d->n_LCC;
        this->n_RyR = d->n_RyR;
        this->MOD = d->MOD;
        this->LCC_STO = d->LCC_STO;
        this->RyR_STO = d->RyR_STO;
        this->rRyRSt = d->rRyRSt;
        this->riRyRSt = d->riRyRSt;
        this->iRyRSt = d->iRyRSt;
        this->cLCCSt = d->cLCCSt;
        this->ifLCCSt = d->ifLCCSt;
        this->if2LCCSt = d->if2LCCSt;
        this->clLCCSt = d->clLCCSt;
        this->iflLCCSt = d->iflLCCSt;
        this->if2lLCCSt = d->if2lLCCSt;
        this->save = d->save;
    }

    static void setV(mreal V){
        CaRU::V = V;
    }
    static void setNai(mreal nai){
        CaRU::Nai = nai;
    }
    static void setSave(bool save){
        CaRU::save = save;
    }
    static void setNLCC(int n_LCC){
        CaRU::n_LCC = n_LCC;
    }
    static void setNRyR(int n_RyR){
        CaRU::n_RyR = n_RyR;
    }
    static void setLCCStochastic(bool s_LCC){
        CaRU::LCC_STO = s_LCC;
    }
    static void setRyRStochastic(bool s_RyR){
        CaRU::RyR_STO = s_RyR;
    }
    static void setModelVariation(int mod){
        CaRU::MOD = mod;
    }

    static void setTUnits(int tUnits){
        CaRU::tUnits = tUnits;
    }

    static void setXUnits(int xUnits){
        CaRU::xUnits = xUnits;
    }

    static void setYUnits(int yUnits){
        CaRU::yUnits = yUnits;
    }

    static void setCaiDiffValue(int id, mreal value){
        CaRU::cais.at(CaRU::getYId(id)).at(CaRU::getXId(id)) = value;
    }

    static mreal getCaiDiffValue(int id){
        return CaRU::cais.at(CaRU::getYId(id)).at(CaRU::getXId(id));
    }

    static void setCaSRDiffValue(int id, mreal value){
        CaRU::casrs.at(CaRU::getYId(id)).at(CaRU::getXId(id)) = value;
    }

    static mreal getCaSRDiffValue(int id){
        return CaRU::casrs.at(CaRU::getYId(id)).at(CaRU::getXId(id));
    }

    static void initiateDefaultInitConditionsCaiDiff(){
        CaRU::cais.clear();
        vetor aux;
        for(int c=0; c<xUnits; c++) aux.push_back(-1.);
        for(int c=0; c<yUnits; c++) CaRU::cais.push_back(aux);
    }

    static void initiateDefaultInitConditionsCaSRDiff(){
        CaRU::casrs.clear();
        vetor aux;
        for(int c=0; c<xUnits; c++) aux.push_back(-1.);
        for(int c=0; c<yUnits; c++) CaRU::casrs.push_back(aux);
    }

    static void saveCaiDiffMatrix(string output_path, mreal t){
        std::ostringstream streamObj3;
        streamObj3 << std::fixed;
        streamObj3 << std::setprecision(0);
        streamObj3 << t;        
        
        string file_name = output_path + "/output_cais_t"+streamObj3.str()+".dat";
        ofstream file (file_name);
        if(file.is_open()){
            for(int y=yUnits-1; y>=0; y--){
                for(int x=0; x<xUnits-1;x++) file<<cais.at(y).at(x)<<"\t";
                file<<cais.at(y).at(xUnits-1)<<endl;
            }
            file.close();
        }
        else cout << "Unable to open file './ca_units/output_cais_t*.dat'"<<endl;
    }

    static void saveCaSRDiffMatrix(string output_path, mreal t){
        std::ostringstream streamObj3;
        streamObj3 << std::fixed;
        streamObj3 << std::setprecision(0);
        streamObj3 << t;

        string file_name = output_path + "/output_casrs_t"+streamObj3.str()+".dat";
        ofstream file (file_name);
        if(file.is_open()){
            for(int y=yUnits-1; y>=0; y--){
                for(int x=0; x<xUnits-1;x++) file<<casrs.at(y).at(x)<<"\t";
                file<<casrs.at(y).at(xUnits-1)<<endl;
            }
            file.close();
        }
        else cout << "Unable to open file './ca_units/output_casrs_t*.dat'"<<endl;
    }

    static void generateCaUnitsPlot(string output_path){
        string plot_ap_path = output_path+"/plot.py";
        ofstream file (plot_ap_path);

        file<<"import matplotlib.pyplot as plt"<<endl;
        file<<"plt.rcParams.update({'font.size': 18})"<<endl;
        file<<"import numpy as np"<<endl;
        file<<"import csv"<<endl;
        file<<endl;
        file<<"##########################################"<<endl;
        file<<"lineSimColor = 'red'"<<endl;
        file<<"#lineOriColor = 'black'"<<endl;
        file<<"##########################################"<<endl;
        file<<endl;
        file<<"##########################################"<<endl;
        file<<"var = []        "<<endl;
        file<<"with open('./output_vars_u0.dat', 'r') as file_v:"<<endl;
        file<<"    lines_v = csv.reader (file_v,delimiter='\t')"<<endl;
        file<<"    for row_v in lines_v:"<<endl;
        file<<"        data = []"<<endl;
        file<<"        # print(row_v[0])"<<endl;
        file<<"        for i in range (len(row_v)):"<<endl;
        file<<"            if(i == 0):"<<endl;
        file<<"                data.append(float(row_v[i])/1000.)"<<endl;
        file<<"            else:"<<endl;
        file<<"                data.append(float(row_v[i]))"<<endl;
        file<<"        var.append(data)"<<endl;
        file<<"var = np.array(var)"<<endl;
        file<<endl;
        file<<"cur = []        "<<endl;
        file<<"with open('./output_curs_u0.dat', 'r') as file_v:"<<endl;
        file<<"    lines_v = csv.reader (file_v,delimiter='\t')"<<endl;
        file<<"    for row_v in lines_v:"<<endl;
        file<<"        data = []"<<endl;
        file<<"        # print(row_v[0])"<<endl;
        file<<"        for i in range (len(row_v)):"<<endl;
        file<<"            if(i == 0):"<<endl;
        file<<"                data.append(float(row_v[i])/1000.)"<<endl;
        file<<"            else:"<<endl;
        file<<"                data.append(float(row_v[i]))"<<endl;
        file<<"        cur.append(data)"<<endl;
        file<<"cur = np.array(cur)"<<endl;
        file<<endl;
        file<<"ext = []        "<<endl;
        file<<"with open('./output_extras_u0.dat', 'r') as file_v:"<<endl;
        file<<"    lines_v = csv.reader (file_v,delimiter='\t')"<<endl;
        file<<"    for row_v in lines_v:"<<endl;
        file<<"        data = []"<<endl;
        file<<"        # print(row_v[0])"<<endl;
        file<<"        for i in range (len(row_v)):"<<endl;
        file<<"            if(i == 0):"<<endl;
        file<<"                data.append(float(row_v[i])/1000.)"<<endl;
        file<<"            else:"<<endl;
        file<<"                data.append(float(row_v[i]))"<<endl;
        file<<"        ext.append(data)"<<endl;
        file<<"ext = np.array(ext)"<<endl;
        file<<endl;
        file<<"#tt = []"<<endl;
        file<<"#with open('../data_TT/output_extras.dat', 'r') as file_e:"<<endl;
        file<<"#    with open('../data_TT/output_vars.dat', 'r') as file_v:"<<endl;
        file<<"#        lines_e = csv.reader (file_e,delimiter='\t')"<<endl;
        file<<"#        lines_v = csv.reader (file_v,delimiter='\t')"<<endl;
        file<<"#        for row_e,row_v in zip(lines_e,lines_v):"<<endl;
        file<<"#            data = []"<<endl;
        file<<"#            for i in range (len(row_e)):"<<endl;
        file<<"#                if(i == 0):"<<endl;
        file<<"#                    data.append(float(row_e[i])/1000.)"<<endl;
        file<<"#                elif(i == 1):"<<endl;
        file<<"#                    data.append(float(row_e[i])*1000.)"<<endl;
        file<<"#                else:"<<endl;
        file<<"#                    data.append(float(row_e[i]))"<<endl;
        file<<"#            for i in range (len(row_v)):"<<endl;
        file<<"#                data.append(float(row_v[i]))"<<endl;
        file<<"#            tt.append(data)"<<endl;
        file<<"#tt = np.array(tt)"<<endl;
        file<<"##########################################"<<endl;
        file<<endl;
        file<<"#IRel#####################################"<<endl;
        file<<"plt.figure()"<<endl;
        file<<"plt.plot(cur[:,0],cur[:,3],color=lineSimColor)"<<endl;
        file<<"#plt.plot(tt[:,0],tt[:,7],color=lineOriColor)"<<endl;
        file<<"plt.xlim(0,0.5)"<<endl;
        file<<"# plt.ylim(0,1)"<<endl;
        file<<"plt.grid(ls=':')"<<endl;
        file<<"plt.xlabel('Time ($s$)')"<<endl;
        file<<"plt.ylabel('$I_{rel}$ ($-$)')"<<endl;
        file<<"plt.savefig('irel.pdf',bbox_inches='tight')"<<endl;
        file<<"##########################################"<<endl;
        file<<endl;
        file<<"#ICaL######################################"<<endl;
        file<<"plt.figure()"<<endl;
        file<<"plt.plot(cur[:,0],cur[:,5],color=lineSimColor)"<<endl;
        file<<"#plt.plot(tt[:,0],tt[:,5],color=lineOriColor)"<<endl;
        file<<"plt.xlim(0,0.5)"<<endl;
        file<<"# plt.ylim(0,1)"<<endl;
        file<<"plt.grid(ls=':')"<<endl;
        file<<"plt.xlabel('Time ($s$)')"<<endl;
        file<<"plt.ylabel('$I_{CaL}$ ($-$)')"<<endl;
        file<<"plt.savefig('ical.pdf',bbox_inches='tight')"<<endl;
        file<<"##########################################"<<endl;
        file<<endl;
        file<<"#OpenICaL#################################"<<endl;
        file<<"plt.figure()"<<endl;
        file<<"plt.plot(ext[:,0],ext[:,1],color=lineSimColor)"<<endl;
        file<<"#plt.plot(tt[:,0],tt[:,5],color=lineOriColor)"<<endl;
        file<<"plt.xlim(0,0.5)"<<endl;
        file<<"# plt.ylim(0,1)"<<endl;
        file<<"plt.grid(ls=':')"<<endl;
        file<<"plt.xlabel('Time ($s$)')"<<endl;
        file<<"plt.ylabel('$O_{CaL}$ ($-$)')"<<endl;
        file<<"plt.savefig('oical.pdf',bbox_inches='tight')"<<endl;
        file<<"##########################################"<<endl;
        file<<endl;
        file<<"#OpenRyR##################################"<<endl;
        file<<"plt.figure()"<<endl;
        file<<"plt.plot(ext[:,0],ext[:,2],color=lineSimColor)"<<endl;
        file<<"#plt.plot(tt[:,0],tt[:,5],color=lineOriColor)"<<endl;
        file<<"plt.xlim(0,0.5)"<<endl;
        file<<"# plt.ylim(0,1)"<<endl;
        file<<"plt.grid(ls=':')"<<endl;
        file<<"plt.xlabel('Time ($s$)')"<<endl;
        file<<"plt.ylabel('$O_{RyR}$ ($-$)')"<<endl;
        file<<"plt.savefig('oryr.pdf',bbox_inches='tight')"<<endl;
        file<<"##########################################"<<endl;
        file<<endl;
        file<<"#Cai######################################"<<endl;
        file<<"plt.figure()"<<endl;
        file<<"plt.plot(ext[:,0],ext[:,3],color=lineSimColor)"<<endl;
        file<<"#plt.plot(tt[:,0],tt[:,5],color=lineOriColor)"<<endl;
        file<<"plt.xlim(0,0.5)"<<endl;
        file<<"# plt.ylim(0,1)"<<endl;
        file<<"plt.grid(ls=':')"<<endl;
        file<<"plt.xlabel('Time ($s$)')"<<endl;
        file<<"plt.ylabel('$[Ca]_{i}$ Free ($-$)')"<<endl;
        file<<"plt.savefig('caifree.pdf',bbox_inches='tight')"<<endl;
        file<<"##########################################"<<endl;

        file.close();
    }

    static void applyDiffusions(mreal dt){
        CaRU::applyCaiDiffusion(dt);
        CaRU::applyCaSRDiffusion(dt);
    }

    static void setDCai(mreal value){
        CaRU::DCai = value;
    }

    static void setDCaSR(mreal value){
        CaRU::DCaSR = value;
    }

private:
    int calcBinomial_1(int t, mreal p);
    int calcBinomial_3(int t, mreal p);
    void initiateDefaultInitConditions();
    void initiateInitConditionsFromFile();
    void initiateDefaultMultiliers();
    void initiateMultiliersFromFile();

    int id;
    vetor y = vetor(_CaRU_numODE_);
    vetor dy = vetor(_CaRU_numODE_);
    vetor mult = vetor(_numM_);
    vector<vetor> save_y;
    vector<vetor> save_c;
    vector<vetor> save_e;
    static vector<vetor> cais;
    static vector<vetor> casrs;

    static bool save;
    static int MOD;                                      //Model variation (TT or TM)
    static int tUnits;                                   // ca_units
    static int xUnits;                                   // ca_units
    static int yUnits;                                   // ca_units
    static int n_LCC;                                    // channels   
    static int n_RyR;                                    // channels
    static bool LCC_STO;                                 //
    static bool RyR_STO;                                 //
    static mreal DCai;
    static mreal DCaSR;
    const static string initCondFilePath;                //
    const static string multFilePath;
    
    // //----Cell Variables-----
        static mreal V;                                      // mV
        static mreal Nai;                                    //
    // //-----------------------

    // //----Model Parameters---
        const mreal F               = 96485.3415;  // C/mmol
        const mreal R               = 8314.472;    // J/molK
        const mreal T               = 310.0;       // K
        const mreal Buf_ss          = 0.4;         //
        const mreal K_buf_ss        = 0.00025;     //
        const mreal Vmax_up         = 0.006375;    //
        const mreal K_up            = 0.00025;     //
        const mreal V_rel           = 0.102;       //
        const mreal V_leak          = 0.00036;     //
        const mreal Buf_c           = 0.2;         //
        const mreal K_buf_c         = 0.001;       //
        const mreal Buf_sr          = 10.;         //
        const mreal K_buf_sr        = 0.3;         //
        const mreal V_c             = 0.016404;    //
        const mreal V_ss            = 0.00005468;  //
        const mreal V_sr            = 0.001094;    //
        const mreal V_xfer          = 0.0038;      //
        const mreal Cm              = 0.185;       //
        const mreal inversevssF2    = 1./(2.*V_ss*F);//
        const mreal Ca_o            = 2.;          //
        const mreal inverseVcF2     = 1./(2.*V_c*F);//
        const mreal RTONF           = (R*T)/F;     //
        const mreal K_pCa           = 0.0005;      //
        const mreal Km_Nai          = 87.5;        //
        const mreal Na_o            = 140.;        //
        const mreal Km_Ca           = 1.38;        //
        const mreal K_sat           = 0.1;         //
        const mreal n               = 0.35;        //
        const mreal k4              = 0.005;       //
        const mreal k3              = 0.06*10.0;   //
        const mreal k2_prime        = 0.045;       //
        const mreal k1_prime        = 0.15*10.0;   //
        const mreal EC              = 1.5;         //
        const mreal max_sr          = 2.5;         //
        const mreal min_sr          = 1.;          //

        const mreal g_CaL           = 0.0000398;   //
        const mreal g_bCa           = 0.000592;    //
        const mreal g_pCa           = 0.1238;      //
        const mreal k_NaCa          = 1000.;       //
    // //------------------------

    // //----Model Variables------         
        mreal yoLCC;
        mreal Cass_free,Cai_free,CaSR_free;
        mreal tau_d,tau_f,tau_f2,tau_fCass;
        mreal d_inf,f_inf,f2_inf,fCass_inf;
        int rRyRSt, riRyRSt, oRyRSt, iRyRSt;                               // units
        int oLCCSt, cLCCSt, clLCCSt, ifLCCSt, iflLCCSt, if2LCCSt, if2lLCCSt;// units
        mreal I_up,I_leak,I_rel,I_xfer,I_CaL,I_bCa,I_pCa,I_NaCa;
        mreal r_C_If,r_C_Cl,r_C_If2;    
        mreal r_If_C,r_If_Ifl;    
        mreal r_If2_C,r_If2_If2l;        
        mreal r_Cl_Ifl,r_Cl_C,r_Cl_If2l,r_Cl_O;
        mreal r_Ifl_O,r_Ifl_Cl,r_Ifl_If;    
        mreal r_If2l_O,r_If2l_Cl,r_If2l_If2;    
        mreal r_O_Ifl,r_O_Cl,r_O_If2l;
        mreal r_R_O, r_O_R, r_R_RI, r_RI_R, r_O_I, r_I_O, r_RI_I, r_I_RI;
    //-------------------------

    static int getXId(int id){
        return id%xUnits;
    }

    static int getYId(int id){
        return  (int) id/yUnits;
    }

    static void applyCaiDiffusion(mreal dt){
        vector<vetor> aux;

        for(int y=0; y<yUnits; y++){
            vetor a;
            for(int x=0; x<yUnits; x++) a.push_back(cais.at(y).at(x));
            aux.push_back(a);
        }

        for(int y=0; y<yUnits; y++){
            for(int x=0; x<yUnits; x++){
                mreal uxy = aux.at(y).at(x);
                mreal uxl1y = (x-1 >= 0) ? aux.at(y).at(x-1) : aux.at(y).at(x);
                mreal uxp1y = (x+1 < xUnits) ? aux.at(y).at(x+1) : aux.at(y).at(x);
                mreal uxyl1 = (y-1 >= 0) ? aux.at(y-1).at(x) : aux.at(y).at(x);
                mreal uxyp1 = (y+1 < yUnits) ? aux.at(y+1).at(x) : aux.at(y).at(x);

                mreal newV = uxy + DCai*dt*(uxp1y-2.0*uxy+uxl1y + uxyp1-2.0*uxy+uxyl1);
                cais.at(y).at(x) = newV;
            }    
        }
    }

    static void applyCaSRDiffusion(mreal dt){
        vector<vetor> aux;

        for(int y=0; y<yUnits; y++){
            vetor a;
            for(int x=0; x<yUnits; x++) a.push_back(casrs.at(y).at(x));
            aux.push_back(a);
        }

        for(int y=0; y<yUnits; y++){
            for(int x=0; x<xUnits; x++){
                mreal uxy = aux.at(y).at(x);
                mreal uxl1y = (x-1 >= 0) ? aux.at(y).at(x-1) : aux.at(y).at(x);
                mreal uxp1y = (x+1 < xUnits) ? aux.at(y).at(x+1) : aux.at(y).at(x);
                mreal uxyl1 = (y-1 >= 0) ? aux.at(y-1).at(x) : aux.at(y).at(x);
                mreal uxyp1 = (y+1 < yUnits) ? aux.at(y+1).at(x) : aux.at(y).at(x);

                mreal newV = uxy + DCaSR*dt*(uxp1y-2.0*uxy+uxl1y + uxyp1-2.0*uxy+uxyl1);
                casrs.at(y).at(x) = newV;
            }    
        }
    }
};

#endif