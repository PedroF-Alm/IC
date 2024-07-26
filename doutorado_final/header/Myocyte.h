#ifndef MYOCYTE_H
#define MYOCYTE_H

#include "CaRU.h"
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <random>
#include <fstream>
#include <chrono>
#include <string>
#include <sys/stat.h>
#include <omp.h>

//Myocyte Model variabels indice definitions-
#define	_Myocyte_numODE_	11
#define _V_             0
#define _Nai_           1
#define _Ki_            2
#define _Xr1_           3
#define _Xr2_           4
#define _Xs_            5
#define _m_             6
#define _h_             7
#define _j_             8
#define _s_             9
#define _r_             10
//--------------------------------------------------

class Myocyte{
public:
    Myocyte();
    Myocyte(int xUnits, int yUnits);
    Myocyte(int xUnits, int yUnits, int n_LCC, int n_RyR);
    Myocyte(int xUnits, int yUnits, int n_LCC, int n_RyR, bool s_LCC, bool s_RyR);
    Myocyte(int xUnits, int yUnits, int n_LCC, int n_RyR, bool s_LCC, bool s_RyR, int mod);
    Myocyte(int xUnits, int yUnits, int n_LCC, int n_RyR, bool s_LCC, bool s_RyR, int mod, mreal DCai, mreal DCaSR);
    void setNLCC(int n_LCC);
    void setNRyR(int n_RyR);
    void setLCCStochastic(bool s_LCC);
    void setRyRStochastic(bool s_RyR);
    void setModelVariation(int mod);
    void solve(mreal dt, mreal t0, mreal tF, mreal printRate, mreal saveRate, int num_threads);
    ~Myocyte();

    static void setOutputPath(string path){
        Myocyte::outputFilePath = path;
        Myocyte::outputCaUnitsFilePath = path+"/ca_units/";
    }

    //----Model Definitions--   
    const static int TT = CaRU::TT;
    const static int TM = CaRU::TM;
    //----------------------- 

    //----Model Parameters---
        const mreal Cm              = 0.185;          //
        const mreal F               = 96485.3415;     // C/mmol
        const mreal R               = 8314.472;       // J/molK
        const mreal T               = 310.0;          // K
        const mreal RTONF           = (R*T)/F;        //
        const mreal V_c             = 0.016404;     //
        const mreal inverseVcF      = 1./(V_c*F);     //
        const mreal P_kna           = 0.03;           // 1
        const mreal K_mk            = 1.;             // mM
        const mreal K_mNa           = 40.;            // mM
        const mreal K_o             = 5.4;          // mM
        const mreal Na_o            = 140.;         // mM

        const mreal stim_period     = 1000.;          // ms
        const mreal stim_start      = 10.;            // ms
        const mreal stim_end        = 30.0e3;       // ms
        const mreal stim_duration   = 1.;             // ms
        const mreal stim_amplitude  = -52.;           // A/F
    //-----------------------

private:
    const static string initCondFilePath;                //
    static string outputFilePath;                  //
    static string outputCaUnitsFilePath;           //
    static int xUnits;
    static int yUnits;
    static int tUnits;
    const static mreal finishMsg;
    
    static bool isPathExist(const string &s){
        struct stat buffer;
        return (stat (s.c_str(), &buffer) == 0);
    }

    static void generatePlotFile(){
        string plot_ap_path = outputFilePath+"/plot.py";
        ofstream file (plot_ap_path.c_str());

        file<<"import matplotlib.pyplot as plt"<<endl;
        file<<"plt.rcParams.update({'font.size': 18})"<<endl;
        file<<"import numpy as np"<<endl;
        file<<"import csv"<<endl;
        file<<endl;
        file<<"##########################################"<<endl;
        file<<"lineSimColor = 'red'"<<endl;
        file<<"lineOriColor = 'black'"<<endl;
        file<<"##########################################"<<endl;
        file<<endl;
        file<<"##########################################"<<endl;
        file<<"var = []        "<<endl;
        file<<"with open('./output_vars.dat', 'r') as file_v:"<<endl;
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
        file<<"with open('./output_curs.dat', 'r') as file_v:"<<endl;
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
        file<<"tt = []"<<endl;
        file<<"with open('../data_TT/output_extras.dat', 'r') as file_e:"<<endl;
        file<<"    with open('../data_TT/output_vars.dat', 'r') as file_v:"<<endl;
        file<<"        lines_e = csv.reader (file_e,delimiter='\t')"<<endl;
        file<<"        lines_v = csv.reader (file_v,delimiter='\t')"<<endl;
        file<<"        for row_e,row_v in zip(lines_e,lines_v):"<<endl;
        file<<"            data = []"<<endl;
        file<<"            for i in range (len(row_e)):"<<endl;
        file<<"                if(i == 0):"<<endl;
        file<<"                    data.append(float(row_e[i])/1000.)"<<endl;
        file<<"                elif(i == 1):"<<endl;
        file<<"                    data.append(float(row_e[i])*1000.)"<<endl;
        file<<"                else:"<<endl;
        file<<"                    data.append(float(row_e[i]))"<<endl;
        file<<"            for i in range (len(row_v)):"<<endl;
        file<<"                data.append(float(row_v[i]))"<<endl;
        file<<"            tt.append(data)"<<endl;
        file<<"tt = np.array(tt)"<<endl;
        file<<"##########################################"<<endl;
        file<<endl;
        file<<"#AP#######################################"<<endl;
        file<<"plt.figure()"<<endl;
        file<<"plt.plot(var[:,0],var[:,1],color=lineSimColor)"<<endl;
        file<<"plt.plot(tt[:,0],tt[:,7],color=lineOriColor)"<<endl;
        file<<"plt.xlim(0,0.5)"<<endl;
        file<<"# plt.ylim(0,1)"<<endl;
        file<<"plt.grid(ls=':')"<<endl;
        file<<"plt.xlabel('Time ($s$)')"<<endl;
        file<<"plt.ylabel('$AP$ ($mV$)')"<<endl;
        file<<"plt.savefig('ap.pdf',bbox_inches='tight')"<<endl;
        file<<"##########################################"<<endl;
        file<<endl;
        file<<"#ICaL######################################"<<endl;
        file<<"plt.figure()"<<endl;
        file<<"plt.plot(cur[:,0],cur[:,5],color=lineSimColor)"<<endl;
        file<<"plt.plot(tt[:,0],tt[:,5],color=lineOriColor)"<<endl;
        file<<"plt.xlim(0,0.5)"<<endl;
        file<<"# plt.ylim(0,1)"<<endl;
        file<<"plt.grid(ls=':')"<<endl;
        file<<"plt.xlabel('Time ($s$)')"<<endl;
        file<<"plt.ylabel('$I_{CaL}$ ($-$)')"<<endl;
        file<<"plt.savefig('ical.pdf',bbox_inches='tight')"<<endl;
        file<<"##########################################"<<endl;

        file.close();

        CaRU::generateCaUnitsPlot(outputCaUnitsFilePath);
    }

    static int getRowPos(int i){
        return int(i/xUnits);
    }

    static int getColPos(int i){
        return int(i % xUnits);
    }
    
    void solveAlgEquations(mreal dt, mreal t);
    void calcDerivatives(mreal dt, mreal t);
    mreal getStim(mreal t);
    mreal getTauH(vetor y);
    mreal getTauJ(vetor y);
    mreal getICaL();
    mreal getIbCa();
    mreal getINaCa();
    mreal getIpCa();
    void initiateDefaultInitConditions();
    void initiateInitConditionsFromFile();
    void solveStep(mreal dt, mreal t);
    void saveVariables();
    void saveCurrents();
    void saveExtras();

    int n_LCC, n_RyR;
    mreal t;
    vector<CaRU*> ca_units;
    mreal caru_ical, caru_inaca, caru_ibca, caru_ipca;
    vetor y = vetor(_Myocyte_numODE_);
    vetor dy = vetor(_Myocyte_numODE_);
    vetor save_t;
    vector<vetor> save_y;
    vector<vetor> save_c;
    vector<vetor> save_e;
    bool toSave;
    int world_size;
    vetor localTUnits;
        
    //----Model Parameters---
        const mreal g_K1                  = 5.405;        // nS/pF
        const mreal g_to                  = 0.073;        // nS/pF
        const mreal g_Kr                  = 0.153;        // nS/pF
        const mreal g_Ks                  = 0.392;        // nS/pF
        const mreal P_NaK                 = 2.724;        // pA/pF
        const mreal g_Na                  = 14.838;       // nS/pF
        const mreal g_bna                 = 0.00029;      // nS/pF
        const mreal g_pK                  = 0.0146;       // nS/pF
    //------------------------

    //----Model Variables------
        mreal I_K1, I_to, I_Kr, I_Ks, I_CaL, I_NaK, I_Na, I_bNa, I_NaCa, I_bCa, I_pK, I_pCa, I_Stim;
        mreal xr1_inf,xr2_inf,xs_inf,m_inf,h_inf,j_inf,s_inf,r_inf;
        mreal tau_xr1,tau_xr2,tau_xs,tau_m,tau_h,tau_j,tau_s,tau_r;
    //-------------------------
};

#endif