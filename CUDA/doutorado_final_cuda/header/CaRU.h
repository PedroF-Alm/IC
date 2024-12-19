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
#include <curand_kernel.h>

using namespace std;

//Type definitions--------
typedef double mreal;
typedef vector<mreal> vetor;
#define Ith(v,i) v[i]
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
    mreal ICal;
    mreal IbCa;
    mreal INaCa;
    mreal IpCa;
} CaruCurs;

class CaRU{
public:
    CaRU(int id);
    ~CaRU();
    __device__ void solveStep(mreal dt, mreal V, mreal Nai, bool save);
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
    mreal *curs;
    mreal *extras;
    mreal *y;
    vector<vetor> save_y;
    vector<vetor> save_c;
    vector<vetor> save_e;
    curandState state;

    //----Model Definitions--   
    const static int TT = 0;
    const static int TM = 1;
    //-----------------------             

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

    static void initiateDefaultInitConditionsCaiDiff(){
        // CaRU::cais.clear();
        // vetor aux;
        // for(int c=0; c<xUnits; c++) aux.push_back(-1.);
        // for(int c=0; c<yUnits; c++) CaRU::cais.push_back(aux);
    }

    static void initiateDefaultInitConditionsCaSRDiff(){
        // CaRU::casrs.clear();
        // vetor aux;
        // for(int c=0; c<xUnits; c++) aux.push_back(-1.);
        // for(int c=0; c<yUnits; c++) CaRU::casrs.push_back(aux);
    }

    static void saveCaiDiffMatrix(mreal *cais, string output_path, mreal t){
        std::ostringstream streamObj3;
        streamObj3 << std::fixed;
        streamObj3 << std::setprecision(0);
        streamObj3 << t;        
        
        string file_name = output_path + "/output_cais_t"+streamObj3.str()+".dat";
        ofstream file (file_name);
        if(file.is_open()){
            for(int y=yUnits-1; y>=0; y--){
                for(int x=0; x<xUnits-1;x++) file<<cais[y * yUnits + x]<<"\t";
                file<<cais[y * yUnits + xUnits - 1]<<endl;
            }
            file.close();
        }
        else cout << "Unable to open file './ca_units/output_cais_t*.dat'"<<endl;
    }

    static void saveCaSRDiffMatrix(mreal *casrs, string output_path, mreal t){
        std::ostringstream streamObj3;
        streamObj3 << std::fixed;
        streamObj3 << std::setprecision(0);
        streamObj3 << t;

        string file_name = output_path + "/output_casrs_t"+streamObj3.str()+".dat";
        ofstream file (file_name);
        if(file.is_open()){
            for(int y=yUnits-1; y>=0; y--){
                for(int x=0; x<xUnits-1;x++) file<<casrs[y * yUnits + x]<<"\t";
                file<<casrs[y * yUnits + xUnits - 1]<<endl;
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

    static void applyDiffusions(mreal *cais, mreal *casrs, mreal dt){
        CaRU::applyCaiDiffusion(cais, dt);
        CaRU::applyCaSRDiffusion(casrs, dt);
    }

    static void setDCai(mreal value){
        CaRU::DCai = value;
    }

    static void setDCaSR(mreal value){
        CaRU::DCaSR = value;
    }

private:
    __device__ void solveAlgEquations(mreal dt, mreal V, mreal Nai);
    __device__ void calcDerivatives(mreal dt, mreal V, mreal Nai);
    __device__ mreal solveQuadratic(mreal a, mreal b, mreal c);
    __device__ void updateRyRs(mreal dt);
    __device__ void updateLCCs(mreal dt);
    __device__ int calcBinomial(int t, mreal p);
    int calcBinomial_1(int t, mreal p);
    int calcBinomial_3(int t, mreal p);
    __device__ int calcBinomial_4(int t, mreal p);
    void initiateDefaultInitConditions();
    void initiateInitConditionsFromFile();
    void initiateDefaultMultiliers();
    void initiateMultiliersFromFile();

    int id;
    mreal *dy;
    mreal *mult;
    
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

    static void applyCaiDiffusion(mreal *cais, mreal dt){
        vector<vetor> aux;

        for(int y=0; y<yUnits; y++){
            vetor a;
            for(int x=0; x<yUnits; x++) a.push_back(cais[y * yUnits + x]);
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
                cais[y * yUnits + x] = newV;
            }    
        }
    }

    static void applyCaSRDiffusion(mreal *casrs, mreal dt){
        vector<vetor> aux;

        for(int y=0; y<yUnits; y++){
            vetor a;
            for(int x=0; x<yUnits; x++) a.push_back(casrs[y * yUnits + x]);
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
                casrs[y * yUnits + x] = newV;
            }    
        }
    }
};


#endif
    
__global__ void solveCaRUStep(CaRU **ca_units, int tUnits, mreal dt, mreal V, mreal Nai, bool save, mreal *cais, mreal *casrs, CaruCurs *caru_curs);