#include "../header/CaRU.h"

mreal CaRU::V       = -85.317147;
mreal CaRU::Nai     = 9.397565;
int CaRU::tUnits    = 1;
int CaRU::xUnits    = 1;
int CaRU::yUnits    = 1;
int CaRU::n_LCC     = 4;
int CaRU::n_RyR     = CaRU::n_LCC*5;
bool CaRU::save     = false;
bool CaRU::LCC_STO  = true;
bool CaRU::RyR_STO  = true;
int CaRU::MOD = CaRU::TM;
mreal CaRU::DCai = 0.0;
mreal CaRU::DCaSR = 0.0;
vector<vetor> CaRU::cais = {{-1.}};
vector<vetor> CaRU::casrs = {{-1.}};
const string CaRU::initCondFilePath = "../input/initCond_CaRU.dat";
const string CaRU::multFilePath = "../input/multipliers.dat";

//----------------------------------------------------------------------

__device__ void solveStep(DataBlock *d, mreal dt){
    calcDerivatives(d, dt);
    for(int c=0; c<_CaRU_numODE_; c++) d->y[c] = d->y[c] + dt*d->dy[c];
    if(d->save){
    
        // d->save_y.push_back(d->y);

    //     mreal curs[8] = {I_up, I_leak, I_rel, I_xfer, I_CaL, I_bCa, I_pCa, I_NaCa};    
    //     save_c.push_back(curs);

    //     vetor extras;
    //     if(!LCC_STO) extras.push_back(yoLCC);
    //     else extras.push_back(((mreal)oLCCSt)/((mreal)n_LCC));
    //     if(!RyR_STO) extras.push_back(Ith(y,_O_));
    //     else extras.push_back(((mreal)oRyRSt)/((mreal)n_RyR));
    //     extras.push_back(Cai_free);
    //     extras.push_back(Cass_free);
    //     extras.push_back(CaSR_free);
    //     save_e.push_back(extras);
    }
}

__global__ void solveCaUnitsStep(SolveStepParameters ssp)
{
    int my_begin = (ssp.tUnits / blockDim.x) * blockIdx.x;
    int my_end = (ssp.tUnits / blockDim.x) * (blockIdx.x + 1);

    while (true)
    {
        if (ssp.V <= -0.5e3)
            break;

        if (threadIdx.x == 0)
        {
            *(ssp.device_V)      = ssp.V;
            *(ssp.device_Nai)    = ssp.Nai;
            *(ssp.device_toSave) = ssp.toSave;
        }

        for (int c = my_begin; c < my_end; c++)
            solveStep(&(ssp.device_c_units[c]), ssp.dt);
        if (ssp.toSave)
        {
            
        }
    }

    printf("Block %d -> The End.", blockIdx.x);
}

CaRU::CaRU(int id){
    // srandom(id*time(NULL));
    srandom(id*2);
    this->id = id;
    this->initiateInitConditionsFromFile();
    this->initiateMultiliersFromFile();
}

void CaRU::initiateDefaultInitConditions(){
    
    Ith(y,_C_)    = 1.0;
    Ith(y,_Cl_)   = 0.0;
    Ith(y,_If_)   = 0.0;
    Ith(y,_Ifl_)  = 0.0;
    Ith(y,_If2_)  = 0.0;
    Ith(y,_If2l_) = 0.0;
    yoLCC         = 1.0 - (Ith(y,_C_) + Ith(y,_Cl_) +Ith(y,_If_) + Ith(y,_Ifl_) + Ith(y,_If2_) + Ith(y,_If2l_));
 
    cLCCSt          = (int)(1.0*((mreal)n_LCC));
    clLCCSt         = (int)(0.0*((mreal)n_LCC));
    ifLCCSt         = (int)(0.0*((mreal)n_LCC));
    iflLCCSt        = (int)(0.0*((mreal)n_LCC));
    if2LCCSt        = (int)(0.0*((mreal)n_LCC));
    if2lLCCSt       = (int)(0.0*((mreal)n_LCC));
    oLCCSt          = n_LCC - (cLCCSt + clLCCSt + ifLCCSt + iflLCCSt + if2LCCSt + if2lLCCSt);

    Ith(y,_R_)    = 0.989568;
    Ith(y,_O_)    = 0.000000;   
    Ith(y,_RI_)   = 0.010423;     
    Ith(y,_I_)    = 0.000000;

    rRyRSt          = (int)(0.989568*((mreal)n_RyR));
    riRyRSt         = (int)(0.010423*((mreal)n_RyR));
    iRyRSt          = (int)(0.000000*((mreal)n_RyR));
    oRyRSt          = n_RyR - (rRyRSt + riRyRSt + iRyRSt);

    Ith(y,_d_)       = 0.000033;
    Ith(y,_f_)       = 0.976018;
    Ith(y,_f2_)      = 0.999494;
    Ith(y,_fCass_)   = 0.999975;
    Ith(y,_CaSR_)    = 12.661646;
    Ith(y,_CaSS_)    = 0.115869;
    Ith(y,_Cai_)     = 0.018842;
}

void CaRU::initiateInitConditionsFromFile(){
    ifstream f;
    string line;
    vetor aux = vetor(_CaRU_numODE_);
    f.open (CaRU::initCondFilePath);
    for(int c=0; c<_CaRU_numODE_; c++){
        getline(f,line);
        aux.at(c) = stod(line);
    }
    f.close();   

    
    Ith(y,_C_)    = Ith(aux,_C_);
    Ith(y,_Cl_)   = Ith(aux,_Cl_);
    Ith(y,_If_)   = Ith(aux,_If_);
    Ith(y,_Ifl_)  = Ith(aux,_Ifl_);
    Ith(y,_If2_)  = Ith(aux,_If2_);
    Ith(y,_If2l_) = Ith(aux,_If2l_);
    yoLCC            = 1.0 - (Ith(y,_C_) + Ith(y,_Cl_) + Ith(y,_If_) + Ith(y,_Ifl_) + Ith(y,_If2_) + Ith(y,_If2l_));
 
    cLCCSt          = (int)(Ith(aux,_C_)*((mreal)n_LCC));
    clLCCSt         = (int)(Ith(aux,_Cl_)*((mreal)n_LCC));
    ifLCCSt         = (int)(Ith(aux,_If_)*((mreal)n_LCC));
    iflLCCSt        = (int)(Ith(aux,_Ifl_)*((mreal)n_LCC));
    if2LCCSt        = (int)(Ith(aux,_If2_)*((mreal)n_LCC));
    if2lLCCSt       = (int)(Ith(aux,_If2l_)*((mreal)n_LCC));
    oLCCSt          = n_LCC - (cLCCSt + clLCCSt + ifLCCSt + iflLCCSt + if2LCCSt + if2lLCCSt);

    Ith(y,_R_)    = Ith(aux,_R_); 
    Ith(y,_O_)    = Ith(aux,_O_);  
    Ith(y,_RI_)   = Ith(aux,_RI_);   
    Ith(y,_I_)    = Ith(aux,_I_); 

    rRyRSt          = (int)(Ith(aux,_R_)*((mreal)n_RyR));
    riRyRSt         = (int)(Ith(aux,_RI_)*((mreal)n_RyR));
    iRyRSt          = (int)(Ith(aux,_I_)*((mreal)n_RyR));
    oRyRSt          = n_RyR - (rRyRSt + riRyRSt + iRyRSt);

    Ith(y,_d_)       = Ith(aux,_d_);
    Ith(y,_f_)       = Ith(aux,_f_);    
    Ith(y,_f2_)      = Ith(aux,_f2_);   
    Ith(y,_fCass_)   = Ith(aux,_fCass_);
    Ith(y,_CaSR_)    = Ith(aux,_CaSR_); 
    Ith(y,_CaSS_)    = Ith(aux,_CaSS_); 
    Ith(y,_Cai_)     = Ith(aux,_Cai_);  
}

void CaRU::initiateDefaultMultiliers(){
    Ith(mult,0) = 1.021870017368191;
    Ith(mult,1) = 1.753962438102943;
    Ith(mult,2) = 0.7663952991525463;
    Ith(mult,3) = 0.632769056408091;
    Ith(mult,4) = 1.7193533070331097;
    Ith(mult,5) = 0.9601108934999598;
    Ith(mult,6) = 0.9378422746139058;
    Ith(mult,7) = 1.2385638121616434;
    Ith(mult,8) = 1.0435388379650659;
    Ith(mult,9) = 1.7717407044148228;
    Ith(mult,10)= 1.721534974901471;
}

void CaRU::initiateMultiliersFromFile(){
    ifstream f;
    string line;
    f.open (CaRU::multFilePath);

    for(int c=0; c<_numM_; c++){
        getline(f,line);
        Ith(mult,c) = stod(line);
    }
    
    f.close();
}

__device__ void solveAlgEquations(DataBlock *d, mreal dt)
{
    mreal V = d->V;
    mreal Nai = d->Nai;

    d->Cai_free = solveQuadratic(-1.0, (d->y[_Cai_] + d->dy[_Cai_] - d->K_buf_c - d->Buf_c), d->K_buf_c * d->y[_Cai_]);
    d->CaSR_free = solveQuadratic(-1.0, (d->y[_CaSR_] + d->dy[_CaSR_] - d->K_buf_sr - d->Buf_sr), d->K_buf_sr * d->y[_CaSR_]);
    d->Cass_free = solveQuadratic(-1.0, (d->y[_CaSS_] + d->y[_CaSS_] - d->K_buf_ss - d->Buf_ss), d->K_buf_ss * d->y[_CaSS_]);

    mreal E_Ca = 0.5 * d->RTONF * (log((d->Ca_o / d->Cai_free)));

    mreal Af = 1102.5 * exp(-(V + 27.) * (V + 27.) / 225.);
    mreal Bf = 200. / (1. + exp((13. - V) / 10.));
    mreal Cf = (180. / (1. + exp((V + 30.) / 10.))) + 20.;
    mreal Af2 = 600. * exp(-(V + 25.) * (V + 25.) / 170.);
    mreal Bf2 = 31. / (1. + exp((25. - V) / 10.));
    mreal Cf2 = 16. / (1. + exp((V + 30.) / 10.));
    mreal Ad = 1.4 / (1. + exp((-35 - V) / 13.)) + 0.25;
    mreal Bd = 1.4 / (1. + exp((V + 5.) / 5.));
    mreal Cd = 1. / (1. + exp((50 - V) / 20.));

    d->tau_d = Ad * Bd + Cd;
    d->tau_f = Af + Bf + Cf;
    d->tau_f2 = Af2 + Bf2 + Cf2;
    d->tau_fCass = 80. / (1. + pow((d->Cass_free / 0.05), 2.)) + 2.;

    d->d_inf = 1. / (1. + exp((-8. - V) / 7.5));
    d->f_inf = 1. / (1. + exp((V + 20.) / 7.));
    d->f2_inf = 0.67 / (1. + exp((V + 35) / 7)) + 0.33;
    d->fCass_inf = 0.6 / (1. + pow((d->Cass_free / 0.05), 2.)) + 0.4;

    mreal oLCC, oRyR;
    updateLCCs(d, dt);
    updateRyRs(d, dt);
    if (d->MOD == d->TT)
    {
        oLCC = d->y[_d_] * d->y[_f_] * d->y[_f2_] * d->y[_fCass_];
        oRyR = d->y[_O_];
    }
    else if (d->MOD == d->TM)
    {
        oLCC = (d->LCC_STO) ? (((mreal)d->oLCCSt) / ((mreal)d->n_LCC)) : d->yoLCC;
        oRyR = (d->RyR_STO) ? (((mreal)d->oRyRSt) / ((mreal)d->n_RyR)) : d->y[_O_];
    }

    d->I_up = d->Vmax_up / (1. + pow(d->K_up, 2.) / pow(d->Cai_free, 2.));
    d->I_leak = d->V_leak * (d->CaSR_free - d->Cai_free);
    d->I_rel = d->V_rel * oRyR * (d->CaSR_free - d->Cass_free);
    d->I_xfer = d->V_xfer * (d->Cass_free - d->Cai_free);
    d->I_CaL = d->g_CaL * oLCC * 4. * (V - 15.) * (d->F * d->F / (d->R * d->T)) * (0.25 * exp(2. * (V - 15.) * d->F / (d->R * d->T)) * d->Cass_free - d->Ca_o) / (exp(2. * (V - 15.) * d->F / (d->R * d->T)) - 1.);
    d->I_bCa = d->g_bCa * (V - E_Ca);
    d->I_pCa = d->g_pCa * d->Cai_free / (d->Cai_free + d->K_pCa);
    d->I_NaCa = d->k_NaCa * (1. / (d->Km_Nai * d->Km_Nai * d->Km_Nai + d->Na_o * d->Na_o * d->Na_o)) * (1. / (d->Km_Ca + d->Ca_o)) * (1. / (1 + d->K_sat * exp((d->n - 1.) * V * d->F / (d->R * d->T)))) * (exp(d->n * V * d->F / (d->R * d->T)) * Nai * Nai * Nai * d->Ca_o - exp((d->n - 1.) * V * d->F / (d->R * d->T)) * d->Na_o * d->Na_o * d->Na_o * d->Cai_free * 2.5);
}

__device__ void calcDerivatives(DataBlock *d, mreal dt){
    solveAlgEquations(d, dt);

    d->dy[_d_] = (d->d_inf - d->y[_d_])/(d->tau_d);
    d->dy[_f_]= (d->f_inf - d->y[_f_])/(d->tau_f);
    d->dy[_f2_] = (d->f2_inf - d->y[_f2_])/(d->tau_f2);
    d->dy[_fCass_] = (d->fCass_inf - d->y[_fCass_])/(d->tau_fCass);
    d->dy[_CaSR_] = (d->I_up-d->I_rel-d->I_leak);
    d->dy[_CaSS_] = (-d->I_xfer*(d->V_c/d->V_ss)+d->I_rel*(d->V_sr/d->V_ss)+(-d->I_CaL*d->inversevssF2*d->Cm));
    d->dy[_Cai_] = ((-(d->I_bCa+d->I_pCa-2.*d->I_NaCa)*d->inverseVcF2*d->Cm)-(d->I_up-d->I_leak)*(d->V_sr/d->V_c)+d->I_xfer);

    d->dy[_C_] = d->r_If_C*d->y[_If_] + d->r_Cl_C*d->y[_Cl_] + d->r_If2_C*d->y[_If2_] - (d->r_C_If + d->r_C_Cl + d->r_C_If2)*d->y[_C_];
    d->dy[_Cl_] = d->r_C_Cl*d->y[_C_] + d->r_Ifl_Cl*d->y[_Ifl_] + d->r_O_Cl*d->yoLCC + d->r_If2l_Cl*d->y[_If2l_] - (d->r_Cl_C + d->r_Cl_Ifl + d->r_Cl_O + d->r_Cl_If2l)*d->y[_Cl_];
    d->dy[_If_] = d->r_C_If*d->y[_C_] + d->r_Ifl_If*d->y[_Ifl_] - (d->r_If_C + d->r_If_Ifl)*d->y[_If_];
    d->dy[_Ifl_] = d->r_If_Ifl*d->y[_If_] + d->r_Cl_Ifl*d->y[_Cl_] + d->r_O_Ifl*d->yoLCC - (d->r_Ifl_If + d->r_Ifl_Cl + d->r_Ifl_O)*d->y[_Ifl_];
    d->dy[_If2_] = d->r_C_If2*d->y[_C_] + d->r_If2l_If2*d->y[_If2l_] - (d->r_If2_C + d->r_If2_If2l)*d->y[_If2_];
    d->dy[_If2l_] = d->r_If2_If2l*d->y[_If2_] + d->r_Cl_If2l*d->y[_Cl_] + d->r_O_If2l*d->yoLCC- (d->r_If2l_If2 + d->r_If2l_Cl + d->r_If2l_O)*d->y[_If2l_];

    d->dy[_R_] = d->r_O_R*d->y[_O_]  + d->r_RI_R*d->y[_RI_]- (d->r_R_O  + d->r_R_RI)*d->y[_R_];
    d->dy[_O_] = d->r_R_O*d->y[_R_]  + d->r_I_O*d->y[_I_]  - (d->r_O_R  + d->r_O_I )*d->y[_O_];
    d->dy[_RI_]= d->r_R_RI*d->y[_R_] + d->r_I_RI*d->y[_I_] - (d->r_RI_R + d->r_RI_I)*d->y[_RI_];
    d->dy[_I_] = d->r_O_I*d->y[_O_]  + d->r_RI_I*d->y[_RI_]- (d->r_I_RI + d->r_I_O )*d->y[_I_];
}

__device__ mreal solveQuadratic(mreal a, mreal b, mreal c){
    mreal delta = b*b - 4.0*a*c;
    mreal x1 = (-b + sqrt(delta))/(2.0*a);
    mreal x2 = (-b - sqrt(delta))/(2.0*a);

    if(x1 >= 0.) return x1;
    else return x2;
}

__device__ void updateRyRs(DataBlock *d, mreal dt){
    mreal kcasr = d->max_sr-((d->max_sr-d->min_sr)/(1.+(d->EC/d->CaSR_free)*(d->EC/d->CaSR_free)));
    mreal k1 = d->k1_prime / kcasr;
    mreal k2 = d->k2_prime * kcasr;

    d->r_R_O = k1*pow(d->Cass_free,2.);
    d->r_R_RI = k2*d->Cass_free;
    d->r_O_R = d->k3;
    d->r_O_I = k2*d->Cass_free;
    d->r_RI_R = d->k4;
    d->r_RI_I = k1*pow(d->Cass_free,2.);
    d->r_I_O = d->k4;
    d->r_I_RI = d->k3;

    if(d->RyR_STO){
        int t_R_O;  
        int t_R_RI; 
        int t_O_R;  
        int t_O_I;  
        int t_RI_I;
        int t_RI_R;            
        int t_I_O;             
        int t_I_RI; 

        bool possible = false;
        while(!possible){
            t_R_O   = calcBinomial(d->rRyRSt, d->r_R_O*dt);
            t_R_RI  = calcBinomial(d->rRyRSt, d->r_R_RI*dt);
            int n_R = d->rRyRSt - (t_R_O + t_R_RI);
            
            t_RI_I  = calcBinomial(d->riRyRSt,d->r_RI_I*dt);
            t_RI_R  = calcBinomial(d->riRyRSt,d->r_RI_R*dt);     
            int n_RI= d->riRyRSt - (t_RI_I + t_RI_R);
            
            t_I_O   = calcBinomial(d->iRyRSt,d->r_I_O*dt);             
            t_I_RI  = calcBinomial(d->iRyRSt,d->r_I_RI*dt); 
            int n_I = d->iRyRSt - (t_I_O + t_I_RI);
            
            t_O_R   = calcBinomial(d->oRyRSt,d->r_O_R*dt); 
            t_O_I   = calcBinomial(d->oRyRSt,d->r_O_I*dt);                 
            int n_O = d->oRyRSt - (t_O_R + t_O_I);
    
            if(n_R < 0 || n_RI < 0 || n_I < 0 || n_O < 0) possible = false;
            else possible = true;

            if(!possible) printf("Valores inválidos!!\t\n");
        }

        int deltaR    = t_RI_R + t_O_R  - (t_R_O  + t_R_RI);
        int deltaRI   = t_R_RI + t_I_RI - (t_RI_I + t_RI_R);
        int deltaI    = t_RI_I + t_O_I  - (t_I_O  + t_I_RI);
        int deltaO    = t_R_O  + t_I_O  - (t_O_R  + t_O_I );

        d->rRyRSt    += deltaR;
        d->riRyRSt   += deltaRI;
        d->iRyRSt    += deltaI;
        d->oRyRSt    = d->n_RyR - (d->rRyRSt + d->riRyRSt + d->iRyRSt);
    }
}

__device__ void updateLCCs(DataBlock *d, mreal dt){
    mreal alpha_d = d->d_inf/d->tau_d;
    mreal beta_f = 1.0*(d->f_inf/d->tau_f);
    mreal beta_f2 = 1.0*(d->f2_inf/d->tau_f2);
    mreal beta_d = (1.-d->d_inf)/d->tau_d;
    mreal alpha_f = 1.0*(1.-d->f_inf)/d->tau_f;
    mreal alpha_f2 = 1.0*(1.-d->f2_inf)/d->tau_f2;

    mreal mahajan_Cass_free;
    mahajan_Cass_free = d->Cass_free;

    mreal cbp = (3.0*d->mult[0]);
    mreal mahajan_f = 1.0/(1.0 + pow(cbp/(mahajan_Cass_free),3.0*d->mult[1]));

    d->r_C_If = d->mult[2]*mahajan_f + alpha_f;
    d->r_C_Cl = alpha_d;
    d->r_C_If2 = alpha_f2;
    
    d->r_If_C = beta_f;
    d->r_If_Ifl = alpha_d;
    
    d->r_If2_C = beta_f2;
    d->r_If2_If2l = alpha_d;
    
    d->r_Cl_Ifl = d->mult[3]*mahajan_f + alpha_f;
    d->r_Cl_C = beta_d;
    d->r_Cl_If2l = alpha_f2;
    d->r_Cl_O = d->mult[4]*2.*alpha_d;
    
    d->r_Ifl_O = d->mult[5]*beta_f;
    d->r_Ifl_Cl = beta_f;
    d->r_Ifl_If = beta_d;
    
    d->r_If2l_O = d->mult[6]*beta_f2;
    d->r_If2l_Cl = beta_f2;
    d->r_If2l_If2 = beta_d;
    
    d->r_O_Ifl = d->mult[7]*mahajan_f + d->mult[8]*3.*alpha_f;
    d->r_O_Cl = d->mult[9]*beta_d;
    d->r_O_If2l = d->mult[10]*1.5*alpha_f2;

    if(d->LCC_STO){
        int t_C_If;
        int t_C_Cl; 
        int t_C_If2; 
        
        int t_If_C; 
        int t_If_Ifl; 
        
        int t_If2_C; 
        int t_If2_If2l; 
        
        int t_Cl_Ifl;
        int t_Cl_C;
        int t_Cl_If2l; 
        int t_Cl_O;
        
        int t_Ifl_O; 
        int t_Ifl_Cl; 
        int t_Ifl_If; 
        
        int t_If2l_O; 
        int t_If2l_Cl; 
        int t_If2l_If2; 
        
        int t_O_Ifl;
        int t_O_Cl; 
        int t_O_If2l;

        bool possible = false;
        while(!possible){
            t_C_If      = calcBinomial(d->cLCCSt,dt*d->r_C_If );
            t_C_Cl      = calcBinomial(d->cLCCSt,dt*d->r_C_Cl );
            t_C_If2     = calcBinomial(d->cLCCSt,dt*d->r_C_If2);
            int n_C     = d->cLCCSt - (t_C_If + t_C_Cl + t_C_If2);
            
            t_If_C      = calcBinomial(d->ifLCCSt,dt*d->r_If_C  );
            t_If_Ifl    = calcBinomial(d->ifLCCSt,dt*d->r_If_Ifl);
            int n_If    = d->ifLCCSt - (t_If_C + t_If_Ifl);
            
            t_If2_C     = calcBinomial(d->if2LCCSt,dt*d->r_If2_C   );
            t_If2_If2l  = calcBinomial(d->if2LCCSt,dt*d->r_If2_If2l);
            int n_If2   = d->if2LCCSt - (t_If2_C + t_If2_If2l);
            
            t_Cl_Ifl    = calcBinomial(d->clLCCSt,dt*d->r_Cl_Ifl );
            t_Cl_C      = calcBinomial(d->clLCCSt,dt*d->r_Cl_C   );
            t_Cl_If2l   = calcBinomial(d->clLCCSt,dt*d->r_Cl_If2l);
            t_Cl_O      = calcBinomial(d->clLCCSt,dt*d->r_Cl_O   );
            int n_Cl    = d->clLCCSt - (t_Cl_Ifl + t_Cl_C + t_Cl_If2l + t_Cl_O);
            
            t_Ifl_O     = calcBinomial(d->iflLCCSt,dt*d->r_Ifl_O );
            t_Ifl_Cl    = calcBinomial(d->iflLCCSt,dt*d->r_Ifl_Cl);
            t_Ifl_If    = calcBinomial(d->iflLCCSt,dt*d->r_Ifl_If);
            int n_Ifl   = d->iflLCCSt - (t_Ifl_O + t_Ifl_Cl + t_Ifl_If);
            
            t_If2l_O    = calcBinomial(d->if2lLCCSt,dt*d->r_If2l_O  );
            t_If2l_Cl   = calcBinomial(d->if2lLCCSt,dt*d->r_If2l_Cl );
            t_If2l_If2  = calcBinomial(d->if2lLCCSt,dt*d->r_If2l_If2);
            int n_If2l  = d->if2lLCCSt - (t_If2l_O + t_If2l_Cl + t_If2l_If2);
            
            t_O_Ifl     = calcBinomial(d->oLCCSt,dt*d->r_O_Ifl );
            t_O_Cl      = calcBinomial(d->oLCCSt,dt*d->r_O_Cl  );
            t_O_If2l    = calcBinomial(d->oLCCSt,dt*d->r_O_If2l);
            int n_O     = d->oLCCSt - (t_O_Ifl + t_O_Cl + t_O_If2l);

            if(n_C < 0 || n_If < 0 || n_If2 < 0 || n_Cl < 0 || n_Ifl < 0 || n_If2l < 0 || n_O < 0) possible = false;
            else possible = true;

            if(!possible) printf("Valores inválidos!!\t\n");
        }

        int deltaC      = t_If_C     + t_If2_C    + t_Cl_C             - (t_C_If   + t_C_Cl     + t_C_If2           );
        int deltaIf     = t_C_If     + t_Ifl_If                        - (t_If_C   + t_If_Ifl                       );
        int deltaIf2    = t_C_If2    + t_If2l_If2                      - (t_If2_C  + t_If2_If2l                     );
        int deltaCl     = t_C_Cl     + t_Ifl_Cl   + t_If2l_Cl + t_O_Cl - (t_Cl_Ifl + t_Cl_C     + t_Cl_If2l + t_Cl_O);
        int deltaIfl    = t_If_Ifl   + t_Cl_Ifl   + t_O_Ifl            - (t_Ifl_O  + t_Ifl_Cl   + t_Ifl_If          );
        int deltaIf2l   = t_If2_If2l + t_Cl_If2l  + t_O_If2l           - (t_If2l_O + t_If2l_Cl  + t_If2l_If2        );
        int deltaO      = t_Cl_O     + t_Ifl_O    + t_If2l_O           - (t_O_Ifl  + t_O_Cl     + t_O_If2l          );

        d->cLCCSt          += deltaC;
        d->ifLCCSt         += deltaIf;
        d->if2LCCSt        += deltaIf2;     
        d->clLCCSt         += deltaCl;
        d->iflLCCSt        += deltaIfl;     
        d->if2lLCCSt       += deltaIf2l;
        d->oLCCSt          = d->n_LCC - (d->cLCCSt+d->ifLCCSt+d->if2LCCSt+d->clLCCSt+d->iflLCCSt+d->if2lLCCSt);
    }
    else d->yoLCC = 1.0 - (d->y[_C_] + d->y[_Cl_] + d->y[_If_] + d->y[_Ifl_] + d->y[_If2_] + d->y[_If2l_]);
}

vetor CaRU::getVariables(){
    return y;
}

vetor CaRU::getDerivatives(){
    return dy;
}

mreal CaRU::getICaL(){
    return I_CaL;
}

mreal CaRU::getIbCa(){
    return I_bCa;
}

mreal CaRU::getINaCa(){
    return I_NaCa;
}

mreal CaRU::getIpCa(){
    return I_pCa;
}

mreal CaRU::getCaiFree(){
    return Cai_free;
}

mreal CaRU::getCaSRFree(){
    return CaSR_free;
}

mreal CaRU::getCai(){
    return Ith(y,_Cai_);
}

mreal CaRU::getCaSR(){
    return Ith(y,_CaSR_);
}

void CaRU::setCai(mreal value){
    Ith(y,_Cai_) = value;
}

void CaRU::setCaSR(mreal value){
    Ith(y,_CaSR_) = value;
}

void CaRU::saveVariables(vetor save_t, string output_path){
    string file_name = output_path + "/output_vars_u"+to_string(id)+".dat";
    ofstream file (file_name);
    if(file.is_open()){
        for(int c=0; c<save_t.size();c++){
            file<<save_t.at(c)<<"\t";
            for(int c1=0; c1<save_y.at(c).size()-1;c1++) file<<save_y.at(c).at(c1)<<"\t";
            file<<save_y.at(c).at(save_y.at(c).size()-1)<<endl;
        }
        file.close();
    }
    else cout << "Unable to open file './ca_units/output_vars_u%d.dat'"<<endl;
}

void CaRU::saveCurrents(vetor save_t, string output_path){
    string file_name = output_path + "/output_curs_u"+to_string(id)+".dat";
    ofstream file (file_name);
    if(file.is_open()){
        for(int c=0; c<save_t.size();c++){
            file<<save_t.at(c)<<"\t";
            for(int c1=0; c1<save_c.at(c).size()-1;c1++) file<<save_c.at(c).at(c1)<<"\t";
            file<<save_c.at(c).at(save_c.at(c).size()-1)<<endl;
        }
        file.close();
    }
    else cout << "Unable to open file './ca_units/output_curs_u%d.dat'"<<endl;
}

void CaRU::saveExtras(vetor save_t, string output_path){
    string file_name = output_path + "/output_extras_u"+to_string(id)+".dat";
    ofstream file (file_name);
    if(file.is_open()){
        for(int c=0; c<save_t.size();c++){
            file<<save_t.at(c)<<"\t";
            for(int c1=0; c1<save_e.at(c).size()-1;c1++) file<<save_e.at(c).at(c1)<<"\t";
            file<<save_e.at(c).at(save_e.at(c).size()-1)<<endl;
        }
        file.close();
    }
    else cout << "Unable to open file './ca_units/output_extras_u%d.dat'"<<endl;
}

//Mudei aqui para a calcBinomial3 pra ficar mais rápido!!!
__device__ int calcBinomial(int t, mreal p){
	// return calcBinomial_1(t,p);
    return calcBinomial_4(t,p);
}
//---------------------------------------------------------

int CaRU::calcBinomial_1(int t, mreal p){
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    default_random_engine generator(seed);
    binomial_distribution<int> bd(t,p);
    return bd(generator);
}

int CaRU::calcBinomial_3(int t, mreal p){
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
	default_random_engine generator(seed);
    uniform_real_distribution<double> distribution(0.0,1.0);
    int res = 0;
    for(int c=0; c<t; c++){
        double prob = distribution(generator);
        if(prob <= p) res++;
    }
    return res;
}

__device__ int calcBinomial_4(int t, mreal p){
    int res = 0;
    for(int c=0; c<t; c++){
        double prob = (double)/*rand()*/ 1/RAND_MAX;
        if(prob <= p) res++;
    }
    return res;
}

void CaRU::printId(){
    cout<<"My Id: "<<this->id<<endl;
}