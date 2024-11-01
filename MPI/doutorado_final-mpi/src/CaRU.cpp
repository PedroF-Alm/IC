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

void CaRU::solveStep(mreal dt){
    calcDerivatives(dt);
    for(int c=0; c<_CaRU_numODE_; c++) Ith(y,c) = Ith(y,c) + dt*Ith(dy,c);
    if(save){
        save_y.push_back(y);

        vetor curs;
        curs.push_back(I_up);
        curs.push_back(I_leak);
        curs.push_back(I_rel);
        curs.push_back(I_xfer);
        curs.push_back(I_CaL);
        curs.push_back(I_bCa);
        curs.push_back(I_pCa);
        curs.push_back(I_NaCa);
        save_c.push_back(curs);

        vetor extras;
        if(!LCC_STO) extras.push_back(yoLCC);
        else extras.push_back(((mreal)oLCCSt)/((mreal)n_LCC));
        if(!RyR_STO) extras.push_back(Ith(y,_O_));
        else extras.push_back(((mreal)oRyRSt)/((mreal)n_RyR));
        extras.push_back(Cai_free);
        extras.push_back(Cass_free);
        extras.push_back(CaSR_free);
        save_e.push_back(extras);
    }
}

void CaRU::solveAlgEquations(mreal dt){
    mreal V = CaRU::V;
    mreal Nai = CaRU::Nai;

    Cai_free = solveQuadratic(-1.0,(Ith(y,_Cai_)+Ith(dy,_Cai_)-K_buf_c-Buf_c),K_buf_c*Ith(y,_Cai_));
    CaSR_free = solveQuadratic(-1.0,(Ith(y,_CaSR_)+Ith(dy,_CaSR_)-K_buf_sr-Buf_sr),K_buf_sr*Ith(y,_CaSR_));
    Cass_free = solveQuadratic(-1.0,(Ith(y,_CaSS_)+Ith(y,_CaSS_)-K_buf_ss-Buf_ss),K_buf_ss*Ith(y,_CaSS_));

    mreal E_Ca = 0.5*RTONF*(log((Ca_o/Cai_free)));

    mreal Af = 1102.5*exp(-(V+27.)*(V+27.)/225.);
    mreal Bf = 200./(1.+exp((13.-V)/10.));
    mreal Cf = (180./(1.+exp((V+30.)/10.)))+20.;
    mreal Af2 = 600.*exp(-(V+25.)*(V+25.)/170.);
    mreal Bf2 = 31./(1.+exp((25.-V)/10.));
    mreal Cf2 = 16./(1.+exp((V+30.)/10.));
    mreal Ad = 1.4/(1.+exp((-35-V)/13.))+0.25;
    mreal Bd = 1.4/(1.+exp((V+5.)/5.));
    mreal Cd = 1./(1.+exp((50-V)/20.));

    tau_d = Ad*Bd + Cd;
    tau_f = Af + Bf + Cf;
    tau_f2 = Af2 + Bf2 + Cf2;
    tau_fCass = 80. / (1. + pow((Cass_free / 0.05),2.)) + 2.;

    d_inf = 1./(1. + exp((-8. - V) / 7.5));
    f_inf = 1./(1. + exp((V + 20.) / 7.));
    f2_inf = 0.67/(1.+exp((V+35)/7))+0.33;
    fCass_inf = 0.6/(1. + pow((Cass_free/0.05),2.)) + 0.4;

    mreal oLCC,oRyR;
    updateLCCs(dt);
    updateRyRs(dt);
    if(MOD == CaRU::TT){
        oLCC = Ith(y,_d_)*Ith(y,_f_)*Ith(y,_f2_)*Ith(y,_fCass_);
        oRyR = Ith(y,_O_);
    }
    else if(MOD == CaRU::TM){
        oLCC = (LCC_STO) ? (((mreal)oLCCSt)/((mreal)n_LCC)) : yoLCC;
        oRyR = (RyR_STO) ? (((mreal)oRyRSt)/((mreal)n_RyR)) : Ith(y,_O_);
    }

    I_up = Vmax_up/(1.+pow(K_up,2.)/pow(Cai_free,2.));
    I_leak = V_leak*(CaSR_free - Cai_free);
    I_rel = V_rel*oRyR*(CaSR_free - Cass_free);
    I_xfer = V_xfer*(Cass_free - Cai_free);
    I_CaL = g_CaL*oLCC*4.*(V-15.)*(F*F/(R*T))*(0.25*exp(2.*(V-15.)*F/(R*T))*Cass_free-Ca_o)/(exp(2.*(V-15.)*F/(R*T))-1.);
    I_bCa = g_bCa*(V-E_Ca);
    I_pCa = g_pCa*Cai_free/(Cai_free+K_pCa);
    I_NaCa = k_NaCa*(1./(Km_Nai*Km_Nai*Km_Nai+Na_o*Na_o*Na_o))*(1./(Km_Ca+Ca_o))*(1./(1+K_sat*exp((n-1.)*V*F/(R*T))))*(exp(n*V*F/(R*T))*Nai*Nai*Nai*Ca_o-exp((n-1.)*V*F/(R*T))*Na_o*Na_o*Na_o*Cai_free*2.5);
}

void CaRU::calcDerivatives(mreal dt){
    solveAlgEquations(dt);

    Ith(dy,_d_) = (d_inf - Ith(y,_d_))/(tau_d);
    Ith(dy,_f_) = (f_inf - Ith(y,_f_))/(tau_f);
    Ith(dy,_f2_) = (f2_inf - Ith(y,_f2_))/(tau_f2);
    Ith(dy,_fCass_) = (fCass_inf - Ith(y,_fCass_))/(tau_fCass);
    Ith(dy,_CaSR_) = (I_up-I_rel-I_leak);
    Ith(dy,_CaSS_) = (-I_xfer*(V_c/V_ss)+I_rel*(V_sr/V_ss)+(-I_CaL*inversevssF2*Cm));
    Ith(dy,_Cai_) = ((-(I_bCa+I_pCa-2.*I_NaCa)*inverseVcF2*Cm)-(I_up-I_leak)*(V_sr/V_c)+I_xfer);

    Ith(dy,_C_) = r_If_C*Ith(y,_If_) + r_Cl_C*Ith(y,_Cl_) + r_If2_C*Ith(y,_If2_) - (r_C_If + r_C_Cl + r_C_If2)*Ith(y,_C_);
    Ith(dy,_Cl_) = r_C_Cl*Ith(y,_C_) + r_Ifl_Cl*Ith(y,_Ifl_) + r_O_Cl*yoLCC + r_If2l_Cl*Ith(y,_If2l_) - (r_Cl_C + r_Cl_Ifl + r_Cl_O + r_Cl_If2l)*Ith(y,_Cl_);
    Ith(dy,_If_) = r_C_If*Ith(y,_C_) + r_Ifl_If*Ith(y,_Ifl_) - (r_If_C + r_If_Ifl)*Ith(y,_If_);
    Ith(dy,_Ifl_) = r_If_Ifl*Ith(y,_If_) + r_Cl_Ifl*Ith(y,_Cl_) + r_O_Ifl*yoLCC - (r_Ifl_If + r_Ifl_Cl + r_Ifl_O)*Ith(y,_Ifl_);
    Ith(dy,_If2_) = r_C_If2*Ith(y,_C_) + r_If2l_If2*Ith(y,_If2l_) - (r_If2_C + r_If2_If2l)*Ith(y,_If2_);
    Ith(dy,_If2l_) = r_If2_If2l*Ith(y,_If2_) + r_Cl_If2l*Ith(y,_Cl_) + r_O_If2l*yoLCC- (r_If2l_If2 + r_If2l_Cl + r_If2l_O)*Ith(y,_If2l_);

    Ith(dy,_R_) = r_O_R*Ith(y,_O_)  + r_RI_R*Ith(y,_RI_)- (r_R_O  + r_R_RI)*Ith(y,_R_);
    Ith(dy,_O_) = r_R_O*Ith(y,_R_)  + r_I_O*Ith(y,_I_)  - (r_O_R  + r_O_I )*Ith(y,_O_);
    Ith(dy,_RI_)= r_R_RI*Ith(y,_R_) + r_I_RI*Ith(y,_I_) - (r_RI_R + r_RI_I)*Ith(y,_RI_);
    Ith(dy,_I_) = r_O_I*Ith(y,_O_)  + r_RI_I*Ith(y,_RI_)- (r_I_RI + r_I_O )*Ith(y,_I_);
}

mreal CaRU::solveQuadratic(mreal a, mreal b, mreal c){
    mreal delta = b*b - 4.0*a*c;
    mreal x1 = (-b + sqrt(delta))/(2.0*a);
    mreal x2 = (-b - sqrt(delta))/(2.0*a);

    if(x1 >= 0.) return x1;
    else return x2;
}

void CaRU::updateRyRs(mreal dt){
    mreal kcasr = max_sr-((max_sr-min_sr)/(1.+(EC/CaSR_free)*(EC/CaSR_free)));
    mreal k1 = k1_prime / kcasr;
    mreal k2 = k2_prime * kcasr;

    r_R_O = k1*pow(Cass_free,2.);
    r_R_RI = k2*Cass_free;
    r_O_R = k3;
    r_O_I = k2*Cass_free;
    r_RI_R = k4;
    r_RI_I = k1*pow(Cass_free,2.);
    r_I_O = k4;
    r_I_RI = k3;

    if(RyR_STO){
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
            t_R_O   = calcBinomial(rRyRSt,r_R_O*dt);
            t_R_RI  = calcBinomial(rRyRSt,r_R_RI*dt);
            int n_R = rRyRSt - (t_R_O + t_R_RI);
            
            t_RI_I  = calcBinomial(riRyRSt,r_RI_I*dt);
            t_RI_R  = calcBinomial(riRyRSt,r_RI_R*dt);     
            int n_RI= riRyRSt - (t_RI_I + t_RI_R);
            
            t_I_O   = calcBinomial(iRyRSt,r_I_O*dt);             
            t_I_RI  = calcBinomial(iRyRSt,r_I_RI*dt); 
            int n_I = iRyRSt - (t_I_O + t_I_RI);
            
            t_O_R   = calcBinomial(oRyRSt,r_O_R*dt); 
            t_O_I   = calcBinomial(oRyRSt,r_O_I*dt);                 
            int n_O = oRyRSt - (t_O_R + t_O_I);
    
            if(n_R < 0 || n_RI < 0 || n_I < 0 || n_O < 0) possible = false;
            else possible = true;

            if(!possible) std::cout<<"Valores inválidos!!\t"<<endl;
        }

        int deltaR    = t_RI_R + t_O_R  - (t_R_O  + t_R_RI);
        int deltaRI   = t_R_RI + t_I_RI - (t_RI_I + t_RI_R);
        int deltaI    = t_RI_I + t_O_I  - (t_I_O  + t_I_RI);
        int deltaO    = t_R_O  + t_I_O  - (t_O_R  + t_O_I );

        rRyRSt    += deltaR;
        riRyRSt   += deltaRI;
        iRyRSt    += deltaI;
        oRyRSt    = n_RyR - (rRyRSt + riRyRSt + iRyRSt);
    }
}

void CaRU::updateLCCs(mreal dt){
    mreal alpha_d = d_inf/tau_d;
    mreal beta_f = 1.0*(f_inf/tau_f);
    mreal beta_f2 = 1.0*(f2_inf/tau_f2);
    mreal beta_d = (1.-d_inf)/tau_d;
    mreal alpha_f = 1.0*(1.-f_inf)/tau_f;
    mreal alpha_f2 = 1.0*(1.-f2_inf)/tau_f2;

    mreal mahajan_Cass_free;
    mahajan_Cass_free = Cass_free;

    mreal cbp = (3.0*Ith(mult,0));
    mreal mahajan_f = 1.0/(1.0 + pow(cbp/(mahajan_Cass_free),3.0*Ith(mult,1)));

    r_C_If = Ith(mult,2)*mahajan_f + alpha_f;
    r_C_Cl = alpha_d;
    r_C_If2 = alpha_f2;
    
    r_If_C = beta_f;
    r_If_Ifl = alpha_d;
    
    r_If2_C = beta_f2;
    r_If2_If2l = alpha_d;
    
    r_Cl_Ifl = Ith(mult,3)*mahajan_f + alpha_f;
    r_Cl_C = beta_d;
    r_Cl_If2l = alpha_f2;
    r_Cl_O = Ith(mult,4)*2.*alpha_d;
    
    r_Ifl_O = Ith(mult,5)*beta_f;
    r_Ifl_Cl = beta_f;
    r_Ifl_If = beta_d;
    
    r_If2l_O = Ith(mult,6)*beta_f2;
    r_If2l_Cl = beta_f2;
    r_If2l_If2 = beta_d;
    
    r_O_Ifl = Ith(mult,7)*mahajan_f + Ith(mult,8)*3.*alpha_f;
    r_O_Cl = Ith(mult,9)*beta_d;
    r_O_If2l = Ith(mult,10)*1.5*alpha_f2;

    if(LCC_STO){
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
            t_C_If      = calcBinomial(cLCCSt,dt*r_C_If );
            t_C_Cl      = calcBinomial(cLCCSt,dt*r_C_Cl );
            t_C_If2     = calcBinomial(cLCCSt,dt*r_C_If2);
            int n_C     = cLCCSt - (t_C_If + t_C_Cl + t_C_If2);
            
            t_If_C      = calcBinomial(ifLCCSt,dt*r_If_C  );
            t_If_Ifl    = calcBinomial(ifLCCSt,dt*r_If_Ifl);
            int n_If    = ifLCCSt - (t_If_C + t_If_Ifl);
            
            t_If2_C     = calcBinomial(if2LCCSt,dt*r_If2_C   );
            t_If2_If2l  = calcBinomial(if2LCCSt,dt*r_If2_If2l);
            int n_If2   = if2LCCSt - (t_If2_C + t_If2_If2l);
            
            t_Cl_Ifl    = calcBinomial(clLCCSt,dt*r_Cl_Ifl );
            t_Cl_C      = calcBinomial(clLCCSt,dt*r_Cl_C   );
            t_Cl_If2l   = calcBinomial(clLCCSt,dt*r_Cl_If2l);
            t_Cl_O      = calcBinomial(clLCCSt,dt*r_Cl_O   );
            int n_Cl    = clLCCSt - (t_Cl_Ifl + t_Cl_C + t_Cl_If2l + t_Cl_O);
            
            t_Ifl_O     = calcBinomial(iflLCCSt,dt*r_Ifl_O );
            t_Ifl_Cl    = calcBinomial(iflLCCSt,dt*r_Ifl_Cl);
            t_Ifl_If    = calcBinomial(iflLCCSt,dt*r_Ifl_If);
            int n_Ifl   = iflLCCSt - (t_Ifl_O + t_Ifl_Cl + t_Ifl_If);
            
            t_If2l_O    = calcBinomial(if2lLCCSt,dt*r_If2l_O  );
            t_If2l_Cl   = calcBinomial(if2lLCCSt,dt*r_If2l_Cl );
            t_If2l_If2  = calcBinomial(if2lLCCSt,dt*r_If2l_If2);
            int n_If2l  = if2lLCCSt - (t_If2l_O + t_If2l_Cl + t_If2l_If2);
            
            t_O_Ifl     = calcBinomial(oLCCSt,dt*r_O_Ifl );
            t_O_Cl      = calcBinomial(oLCCSt,dt*r_O_Cl  );
            t_O_If2l    = calcBinomial(oLCCSt,dt*r_O_If2l);
            int n_O     = oLCCSt - (t_O_Ifl + t_O_Cl + t_O_If2l);

            if(n_C < 0 || n_If < 0 || n_If2 < 0 || n_Cl < 0 || n_Ifl < 0 || n_If2l < 0 || n_O < 0) possible = false;
            else possible = true;

            if(!possible) std::cout<<"Valores inválidos!!\t"<<endl;
        }

        int deltaC      = t_If_C     + t_If2_C    + t_Cl_C             - (t_C_If   + t_C_Cl     + t_C_If2           );
        int deltaIf     = t_C_If     + t_Ifl_If                        - (t_If_C   + t_If_Ifl                       );
        int deltaIf2    = t_C_If2    + t_If2l_If2                      - (t_If2_C  + t_If2_If2l                     );
        int deltaCl     = t_C_Cl     + t_Ifl_Cl   + t_If2l_Cl + t_O_Cl - (t_Cl_Ifl + t_Cl_C     + t_Cl_If2l + t_Cl_O);
        int deltaIfl    = t_If_Ifl   + t_Cl_Ifl   + t_O_Ifl            - (t_Ifl_O  + t_Ifl_Cl   + t_Ifl_If          );
        int deltaIf2l   = t_If2_If2l + t_Cl_If2l  + t_O_If2l           - (t_If2l_O + t_If2l_Cl  + t_If2l_If2        );
        int deltaO      = t_Cl_O     + t_Ifl_O    + t_If2l_O           - (t_O_Ifl  + t_O_Cl     + t_O_If2l          );

        cLCCSt          += deltaC;
        ifLCCSt         += deltaIf;
        if2LCCSt        += deltaIf2;     
        clLCCSt         += deltaCl;
        iflLCCSt        += deltaIfl;     
        if2lLCCSt       += deltaIf2l;
        oLCCSt          = n_LCC - (cLCCSt+ifLCCSt+if2LCCSt+clLCCSt+iflLCCSt+if2lLCCSt);
    }
    else yoLCC = 1.0 - (Ith(y,_C_) + Ith(y,_Cl_) + Ith(y,_If_) + Ith(y,_Ifl_) + Ith(y,_If2_) + Ith(y,_If2l_));
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
int CaRU::calcBinomial(int t, mreal p){
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

int CaRU::calcBinomial_4(int t, mreal p){
    int res = 0;
    for(int c=0; c<t; c++){
        double prob = (double)rand()/RAND_MAX;
        if(prob <= p) res++;
    }
    return res;
}

void CaRU::printId(){
    cout<<"My Id: "<<this->id<<endl;
}