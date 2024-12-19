#include "Myocyte.h"

const string Myocyte::initCondFilePath = "../input/initCond_Myocyte.dat";
string Myocyte::outputFilePath = "../exit/lastSim/";
string Myocyte::outputCaUnitsFilePath = outputFilePath+"/ca_units/";
int Myocyte::xUnits = 1;
int Myocyte::yUnits = 1;
int Myocyte::tUnits = 1;
const mreal Myocyte::finishMsg = -1.0e3;
//------------------------------------------------------------------

Myocyte::Myocyte(){
    this->initiateInitConditionsFromFile();

    CaRU::setTUnits(tUnits);
    CaRU::setXUnits(xUnits);
    CaRU::setYUnits(yUnits);
    CaRU::initiateDefaultInitConditionsCaiDiff();
    CaRU::initiateDefaultInitConditionsCaSRDiff();
    //for(int c=0; c<tUnits; c++) this->ca_units.push_back(new CaRU(c));
}

Myocyte::Myocyte(int xUnits, int yUnits){
    this->initiateInitConditionsFromFile();
    Myocyte::xUnits = xUnits;
    Myocyte::yUnits = yUnits;
    Myocyte::tUnits = (xUnits*yUnits);

    CaRU::setTUnits(tUnits);
    CaRU::setXUnits(xUnits);
    CaRU::setYUnits(yUnits);
    CaRU::initiateDefaultInitConditionsCaiDiff();
    CaRU::initiateDefaultInitConditionsCaSRDiff();
    //for(int c=0; c<tUnits; c++) this->ca_units.push_back(new CaRU(c));
}

Myocyte::Myocyte(int xUnits, int yUnits, int n_LCC, int n_RyR){
    this->initiateInitConditionsFromFile();
    Myocyte::xUnits = xUnits;
    Myocyte::yUnits = yUnits;
    this->n_LCC = n_LCC;
    this->n_RyR = n_RyR;
    Myocyte::tUnits = (xUnits*yUnits);

    CaRU::setTUnits(tUnits);
    CaRU::setXUnits(xUnits);
    CaRU::setYUnits(yUnits);
    CaRU::initiateDefaultInitConditionsCaiDiff();
    CaRU::initiateDefaultInitConditionsCaSRDiff();
    setNLCC(n_LCC);
    setNRyR(n_RyR);
    //for(int c=0; c<tUnits; c++) this->ca_units.push_back(new CaRU(c));
}

Myocyte::Myocyte(int xUnits, int yUnits, int n_LCC, int n_RyR, bool s_LCC, bool s_RyR){
    this->initiateInitConditionsFromFile();
    Myocyte::xUnits = xUnits;
    Myocyte::yUnits = yUnits;
    this->n_LCC = n_LCC;
    this->n_RyR = n_RyR;
    Myocyte::tUnits = (xUnits*yUnits);

    CaRU::setTUnits(tUnits);
    CaRU::setXUnits(xUnits);
    CaRU::setYUnits(yUnits);
    CaRU::initiateDefaultInitConditionsCaiDiff();
    CaRU::initiateDefaultInitConditionsCaSRDiff();
    setNLCC(n_LCC);
    setNRyR(n_RyR);
    setLCCStochastic(s_LCC);
    setRyRStochastic(s_RyR);
    //for(int c=0; c<tUnits; c++) this->ca_units.push_back(new CaRU(c));
}

Myocyte::Myocyte(int xUnits, int yUnits, int n_LCC, int n_RyR, bool s_LCC, bool s_RyR, int mod){
    this->initiateInitConditionsFromFile();
    Myocyte::xUnits = xUnits;
    Myocyte::yUnits = yUnits;
    this->n_LCC = n_LCC;
    this->n_RyR = n_RyR;
    Myocyte::tUnits = (xUnits*yUnits);

    CaRU::setTUnits(tUnits);
    CaRU::setXUnits(xUnits);
    CaRU::setYUnits(yUnits);
    CaRU::initiateDefaultInitConditionsCaiDiff();
    CaRU::initiateDefaultInitConditionsCaSRDiff();
    setNLCC(n_LCC);
    setNRyR(n_RyR);
    setLCCStochastic(s_LCC);
    setRyRStochastic(s_RyR);
    setModelVariation(mod);
    //for(int c=0; c<tUnits; c++) this->ca_units.push_back(new CaRU(c));
}

Myocyte::Myocyte(int xUnits, int yUnits, int n_LCC, int n_RyR, bool s_LCC, bool s_RyR, int mod, mreal DCai, mreal DCaSR){
    this->initiateInitConditionsFromFile();
    Myocyte::xUnits = xUnits;
    Myocyte::yUnits = yUnits;
    this->n_LCC = n_LCC;
    this->n_RyR = n_RyR;
    Myocyte::tUnits = (xUnits*yUnits);

    CaRU::setTUnits(tUnits);
    CaRU::setXUnits(xUnits);
    CaRU::setYUnits(yUnits);
    CaRU::initiateDefaultInitConditionsCaiDiff();
    CaRU::initiateDefaultInitConditionsCaSRDiff();
    CaRU::setDCai(DCai);
    CaRU::setDCaSR(DCaSR);
    setNLCC(n_LCC);
    setNRyR(n_RyR);
    setLCCStochastic(s_LCC);
    setRyRStochastic(s_RyR);
    setModelVariation(mod);
    //for(int c=0; c<tUnits; c++) this->ca_units.push_back(new CaRU(c));
}

void Myocyte::initiateDefaultInitConditions(){
    Ith(y,_V_)    = -85.317147;
    Ith(y,_Nai_)  = 9.397565;
    Ith(y,_Ki_)   = 135.388663;
    Ith(y,_Xr1_)  = 0.000210;
    Ith(y,_Xr2_)  = 0.472082;
    Ith(y,_Xs_)   = 0.003250;
    Ith(y,_m_)    = 0.001684;
    Ith(y,_h_)    = 0.747286;
    Ith(y,_j_)    = 0.746880;
    Ith(y,_s_)    = 0.633555;
    Ith(y,_r_)    = 0.000000;
}

void Myocyte::initiateInitConditionsFromFile(){
    ifstream f;
    string line;
    vetor aux = vetor(_CaRU_numODE_);
    f.open (Myocyte::initCondFilePath.c_str());
    for(int c=0; c<_Myocyte_numODE_; c++){
        getline(f,line);
        aux.at(c) = stod(line);
    }
    f.close(); 

    Ith(y,_V_)    = Ith(aux,_V_);
    Ith(y,_Nai_)  = Ith(aux,_Nai_);
    Ith(y,_Ki_)   = Ith(aux,_Ki_);
    Ith(y,_Xr1_)  = Ith(aux,_Xr1_);
    Ith(y,_Xr2_)  = Ith(aux,_Xr2_);
    Ith(y,_Xs_)   = Ith(aux,_Xs_);
    Ith(y,_m_)    = Ith(aux,_m_);
    Ith(y,_h_)    = Ith(aux,_h_);
    Ith(y,_j_)    = Ith(aux,_j_);
    Ith(y,_s_)    = Ith(aux,_s_);
    Ith(y,_r_)    = Ith(aux,_r_);
}

void Myocyte::solveStep(mreal dt, mreal t){
    this->calcDerivatives(dt,t);
    for(int c=0; c<_Myocyte_numODE_; c++) Ith(y,c) = Ith(y,c) + dt*Ith(dy,c);
}

void Myocyte::solve(mreal dt, mreal t0, mreal tF, mreal printRate, mreal saveRate){
    cout << "Creating output dir." << endl;
    string command = "rm -r ";
    if (isPathExist(outputFilePath)){
        command += outputFilePath;
        system(command.c_str());
    }
    command = "mkdir " + outputFilePath;
    system(command.c_str());
    cout << "Output dir created." << endl;
    command = "mkdir " + outputCaUnitsFilePath;
    system(command.c_str());
    cout << "CaRU units dir created." << endl;
    cout << "Beginning solve Myocyte model." << endl;
    t = t0;

    caru_ical = 0.0;
    caru_inaca = 0.0;
    caru_ibca = 0.0;
    caru_ipca = 0.0;

    for (int c = 0; c < tUnits; c++) ca_units.push_back(new CaRU(c));

    double caru_solve_time = 0, caru_diff_time = 0, myocyte_solve_time = 0;
    auto begin = chrono::high_resolution_clock::now(), end = chrono::high_resolution_clock::now();      
    double pBegin, pEnd;

    for (int c = 0; t < tF; c++)
    {
        toSave = (saveRate > 0.0 and c % (int)((1.0 / saveRate)) == 0);
        if (printRate > 0.0 and c % (int)((1.0 / printRate)) == 0)
            cout << "Grid: " << xUnits << "x" << yUnits << "\t n_LCC: " << n_LCC << "\t" << "Time (ms):" << t << endl;
        if (toSave)
        {
            save_t.push_back(t);
            save_y.push_back(y);
        }

        for (int u = 0; u < tUnits; u++)
        {
            CaRU::setCaiDiffValue(u, ca_units.at(u)->getCai());
            CaRU::setCaSRDiffValue(u, ca_units.at(u)->getCaSR());
            if (CaRU::getCaiDiffValue(u) > -1.)
                ca_units.at(u)->setCai(CaRU::getCaiDiffValue(u));
            if (CaRU::getCaSRDiffValue(u) > -1.)
                ca_units.at(u)->setCaSR(CaRU::getCaSRDiffValue(u));
        }

        CaRU::setV(Ith(y, _V_));
        CaRU::setNai(Ith(y, _Nai_));
        CaRU::setSave(toSave);

        pBegin = omp_get_wtime();

        #pragma omp for schedule(static, tUnits / omp_get_num_threads())
        for (int c = 0; c < tUnits; c++) {
            ca_units.at(c)->solveStep(dt);
        }

        pEnd = omp_get_wtime();

        caru_solve_time += pEnd - pBegin;

        begin = chrono::high_resolution_clock::now();
        solveStep(dt, t);
        end = chrono::high_resolution_clock::now();

        myocyte_solve_time += chrono::duration_cast<chrono::nanoseconds>(end - begin).count() * 1e-9;

        begin = chrono::high_resolution_clock::now();
        CaRU::applyDiffusions(dt);
        end = chrono::high_resolution_clock::now();

        caru_diff_time += chrono::duration_cast<chrono::nanoseconds>(end - begin).count() * 1e-9;

        if (toSave)
        {
            vetor curs;
            curs.push_back(I_K1);
            curs.push_back(I_to);
            curs.push_back(I_Kr);
            curs.push_back(I_Ks);
            curs.push_back(I_CaL);
            curs.push_back(I_NaK);
            curs.push_back(I_Na);
            curs.push_back(I_bNa);
            curs.push_back(I_NaCa);
            curs.push_back(I_bCa);
            curs.push_back(I_pK);
            curs.push_back(I_pCa);
            curs.push_back(I_Stim);
            save_c.push_back(curs);

            CaRU::saveCaiDiffMatrix(outputCaUnitsFilePath, t);
            CaRU::saveCaSRDiffMatrix(outputCaUnitsFilePath, t);
        }

        // if(saveRate > 0.0 and c%(int)((1.0/saveRate)) == 0){
        //     vetor extras;
        //     save_e.push_back(extras);
        // }
        t += dt;
    }

    cout << "Fineshed solving." << endl;

    if (saveRate > 0.0){
        // Myocyte savings
        cout << "Saving variables results." << endl;
        saveVariables();
        cout << "Fineshed saving." << endl;
        cout << "Saving Currents results." << endl;
        saveCurrents();
        cout << "Fineshed saving." << endl;
        cout << "Saving Extra results." << endl;
        saveExtras();
        cout << "Fineshed saving." << endl;
        cout << "Plotting results." << endl;
        generatePlotFile();
        command = "(cd " + outputFilePath + "/;python plot.py)";
        system(command.c_str());
        command = "(cd " + outputCaUnitsFilePath + "/;python plot.py)";
        system(command.c_str());
        cout << "Fineshed plotting." << endl;

        // CaRUs savings
        cout << "Saving variables results." << endl;
        for (int c = 0; c < tUnits; c++) ca_units.at(c)->saveVariables(save_t, outputCaUnitsFilePath);
        cout << "Fineshed saving." << endl;
        cout << "Saving Currents results." << endl;
        for (int c = 0; c < tUnits; c++) ca_units.at(c)->saveCurrents(save_t, outputCaUnitsFilePath);
        cout << "Fineshed saving." << endl;
        cout << "Saving Extra results." << endl;
        for (int c = 0; c < tUnits; c++) ca_units.at(c)->saveExtras(save_t, outputCaUnitsFilePath);
        cout << "Fineshed saving." << endl;
    }

    for (int c = 0; c < tUnits; c++) delete ca_units.at(c);

    ofstream file;
    file.open("/home/lince/Documentos/IC/testes/tempos/openmp.txt", ios::app);
    file.precision(8);
    file << "\nCaRUs solving time: " << fixed << caru_solve_time << "s\n";
    file << "Myocyte solving time: " << fixed << myocyte_solve_time << "s\n";
    file << "CaRUs diffusion time: " << fixed << caru_diff_time << "s\n\n";
    file.close();
}

void Myocyte::solveAlgEquations(mreal dt, mreal t){
    mreal E_K = RTONF*(log((K_o/Ith(y,_Ki_))));
    mreal E_Ks = RTONF*(log((K_o + P_kna*Na_o)/(Ith(y,_Ki_) + P_kna*Ith(y,_Nai_))));
    mreal E_Na = RTONF*(log((Na_o/Ith(y,_Nai_))));

    mreal Ak1 = 0.1/(1.+exp(0.06*(Ith(y,_V_)-E_K-200.)));
    mreal Bk1 = (3.*exp(0.0002*(Ith(y,_V_)-E_K+100.)) + exp(0.1*(Ith(y,_V_)-E_K-10.)))/(1.+exp(-0.5*(Ith(y,_V_)-E_K)));
    mreal rec_iK1 = Ak1/(Ak1+Bk1);
    mreal rec_iNaK = (1./(1.+0.1245*exp(-0.1*Ith(y,_V_)*F/(R*T))+0.0353*exp(-Ith(y,_V_)*F/(R*T))));
    mreal rec_ipK =1./(1.+exp((25. - Ith(y,_V_))/5.98));

    mreal alpha_xr1 = 450./(1.+exp((-45.-Ith(y,_V_))/10.));
    mreal beta_xr1 = 6./(1.+exp((Ith(y,_V_) - (-30.)) / 11.5));
    mreal alpha_xr2 = 3./(1.+exp((-60.-Ith(y,_V_)) / 20.));
    mreal beta_xr2 = 1.12/(1.+exp((Ith(y,_V_)-60.)/20.));
    mreal alpha_xs = (1400./(sqrt(1.+exp((5.-Ith(y,_V_))/6.))));
    mreal beta_xs = (1./(1.+exp((Ith(y,_V_)-35.)/15.)));
    mreal alpha_m = 1./(1.+exp((-60.-Ith(y,_V_))/5.));
    mreal beta_m = 0.1/(1.+exp((Ith(y,_V_)+35.)/5.))+0.1/(1.+exp((Ith(y,_V_)-50.) / 200.));

    tau_xr1 = alpha_xr1*beta_xr1;
    tau_xr2 = alpha_xr2*beta_xr2;
    tau_xs = alpha_xs*beta_xs + 80.;
    tau_m = alpha_m*beta_m;
    tau_h = getTauH(y);
    tau_j = getTauJ(y);    
    tau_s = 1000.*exp(-pow((Ith(y,_V_)+67.),2.)/1000.)+8.;
    tau_r = 9.5*exp(-(Ith(y,_V_)+40.)*(Ith(y,_V_)+40.)/1800.)+0.8;

    xr1_inf = 1./(1.+exp((-26.-Ith(y,_V_))/7.));
    xr2_inf = 1./(1.+exp((Ith(y,_V_) - (-88.)) / 24.));
    xs_inf = 1./(1.+exp((-5.-Ith(y,_V_))/14.));
    m_inf = 1./((1.+exp((-56.86-Ith(y,_V_))/9.03))*(1.+exp((-56.86-Ith(y,_V_))/9.03)));
    h_inf = 1./((1.+exp((Ith(y,_V_)+71.55)/7.43))*(1.+exp((Ith(y,_V_)+71.55)/7.43)));
    j_inf = h_inf;
    s_inf = 1./(1.+exp((Ith(y,_V_)+28.)/5.));   
    r_inf = 1./(1.+exp((20.-Ith(y,_V_))/6.));

    I_K1 = g_K1*rec_iK1*(Ith(y,_V_) - E_K);
    I_to = g_to*Ith(y,_r_)*Ith(y,_s_)*(Ith(y,_V_) - E_K);
    I_Kr = g_Kr*sqrt(K_o/5.4)*Ith(y,_Xr1_)*Ith(y,_Xr2_)*(Ith(y,_V_) - E_K);
    I_Ks = g_Ks*pow(Ith(y,_Xs_),2.)*(Ith(y,_V_) - E_Ks);
    I_CaL = getICaL();
    I_NaK = P_NaK*(K_o/(K_o + K_mk))*(Ith(y,_Nai_)/(Ith(y,_Nai_)+K_mNa))*rec_iNaK;
    I_Na = g_Na*pow(Ith(y,_m_),3.)*Ith(y,_h_)*Ith(y,_j_)*(Ith(y,_V_) - E_Na);
    I_bNa = g_bna*(Ith(y,_V_) - E_Na);
    I_NaCa = getINaCa();
    I_bCa = getIbCa();
    I_pK = g_pK*rec_ipK*(Ith(y,_V_)-E_K);
    I_pCa = getIpCa();
    I_Stim = getStim(t);
}

void Myocyte::calcDerivatives(mreal dt, mreal t){
    // CaRU::setV(Ith(y,_V_));
    // CaRU::setNai(Ith(y,_Nai_));
    //for(int c=0; c<tUnits; c++) ca_units.at(c)->solveStep(dt);

    solveAlgEquations(dt,t);

    Ith(dy,_V_) = -(I_K1 + I_to + I_Kr + I_Ks + I_CaL + I_NaK + I_Na + I_bNa + I_NaCa + I_bCa + I_pK + I_pCa + I_Stim);
    Ith(dy,_Nai_) = -(I_Na+I_bNa+3.*I_NaK+3.*I_NaCa)*inverseVcF*Cm;
    Ith(dy,_Ki_) = -(I_Stim+I_K1+I_to+I_Kr+I_Ks-2.*I_NaK+I_pK)*inverseVcF*Cm;
    Ith(dy,_Xr1_) = (xr1_inf - Ith(y,_Xr1_))/(tau_xr1);
    Ith(dy,_Xr2_) = (xr2_inf - Ith(y,_Xr2_))/(tau_xr2);
    Ith(dy,_Xs_) = (xs_inf - Ith(y,_Xs_))/(tau_xs);
    Ith(dy,_m_) = (m_inf - Ith(y,_m_))/(tau_m);
    Ith(dy,_h_) = (h_inf - Ith(y,_h_))/(tau_h);
    Ith(dy,_j_) = (j_inf - Ith(y,_j_))/(tau_j);
    Ith(dy,_s_) = (s_inf - Ith(y,_s_))/(tau_s);
    Ith(dy,_r_) = (r_inf - Ith(y,_r_))/(tau_r);
}

mreal Myocyte::getICaL(){
    mreal all_Curs = 0.0;

    for(int c=0; c<tUnits; c++) all_Curs += ca_units.at(c)->getICaL();
    all_Curs /= (mreal)tUnits;

    return all_Curs;
}

mreal Myocyte::getIbCa(){
    mreal all_Curs = 0.0;

    for(int c=0; c<tUnits; c++) all_Curs += ca_units.at(c)->getIbCa();
    all_Curs /= (mreal)tUnits;

    return all_Curs;
}

mreal Myocyte::getINaCa(){
    mreal all_Curs = 0.0;

    for(int c=0; c<tUnits; c++) all_Curs += ca_units.at(c)->getINaCa();
    all_Curs /= (mreal)tUnits;

    return all_Curs;
}

mreal Myocyte::getIpCa(){
    mreal all_Curs = 0.0;

    for(int c=0; c<tUnits; c++) all_Curs += ca_units.at(c)->getIpCa();
    all_Curs /= (mreal)tUnits;

    return all_Curs;
}

mreal Myocyte::getStim(mreal t){
    if(((t >= stim_start) && (t <= stim_end) && (((t - stim_start) - (floor(((t - stim_start)/stim_period))*stim_period)) <= stim_duration))) {
        return (stim_amplitude);
    }
    else {
        return (0.0);
    }
}

mreal Myocyte::getTauH(vetor y){
    if (Ith(y,_V_)>=-40.){
        mreal AH_1 = 0.; 
        mreal BH_1 = (0.77/(0.13*(1.+exp(-(Ith(y,_V_)+10.66)/11.1))));
        return 1.0/(AH_1+BH_1);
    }
    else{
        mreal AH_2 = (0.057*exp(-(Ith(y,_V_)+80.)/6.8));
        mreal BH_2 = (2.7*exp(0.079*Ith(y,_V_))+(3.1e5)*exp(0.3485*Ith(y,_V_)));
        return 1.0/(AH_2+BH_2);
    }
}

mreal Myocyte::getTauJ(vetor y){
    if(Ith(y,_V_)>=-40.){
        mreal AJ_1 = 0.;      
        mreal BJ_1 = (0.6*exp((0.057)*Ith(y,_V_))/(1.+exp(-0.1*(Ith(y,_V_)+32.))));
        return  1.0/(AJ_1+BJ_1);
    }
    else{
        mreal AJ_2 = (((-2.5428e4)*exp(0.2444*Ith(y,_V_))-(6.948e-6)*exp(-0.04391*Ith(y,_V_)))*(Ith(y,_V_)+37.78)/(1.+exp(0.311*(Ith(y,_V_)+79.23))));
        mreal BJ_2 = (0.02424*exp(-0.01052*Ith(y,_V_))/(1.+exp(-0.1378*(Ith(y,_V_)+40.14))));
        return 1.0/(AJ_2+BJ_2);
    }
}

void Myocyte::setNLCC(int n_LCC){
    CaRU::setNLCC(n_LCC);
}

void Myocyte::setNRyR(int n_RyR){
    CaRU::setNRyR(n_RyR);
}

void Myocyte::setLCCStochastic(bool s_LCC){
    CaRU::setLCCStochastic(s_LCC);
}

void Myocyte::setRyRStochastic(bool s_RyR){
    CaRU::setRyRStochastic(s_RyR);
}

void Myocyte::setModelVariation(int mod){
    CaRU::setModelVariation(mod);
}

void Myocyte::saveVariables(){
    string vars_path = outputFilePath+"/output_vars.dat";
    ofstream file (vars_path.c_str());
    if(file.is_open()){
        for(int c=0; c<save_t.size();c++){
            file<<save_t.at(c)<<"\t";
            for(int c1=0; c1<save_y.at(c).size()-1;c1++) file<<save_y.at(c).at(c1)<<"\t";
            file<<save_y.at(c).at(save_y.at(c).size()-1)<<endl;
        }
        file.close();
    }
    else cout << "Unable to open file './output_vars.dat'"<<endl;

    /*// char command[100];
    sprintf(command,"mkdir %s/cai_X_time",output_path.c_str());
    system(command);
    for(int t=0; t<save_t.size(); t++){
        char filename[100];
        sprintf(filename,"%s/cai_X_time/output_t%d.dat",output_path.c_str(),t);
        file.clear();file.open(filename);
        if(file.is_open()){
            for(int u=0; u<tUnits;u++){
                file<<(ca_units.at(u)->getCai(t)*1.0e3);
                if((u+1)%xUnits == 0) file<<"\n";
                else file<<"\t";
            }
            file.close();
        }
        else cout << "Unable to open file './exit/cai_X_time/output_uX.dat'"<<endl;
    }*/
}

void Myocyte::saveCurrents(){
    string curs_path = outputFilePath+"/output_curs.dat";
    ofstream file (curs_path.c_str());
    if(file.is_open()){
        for(int c=0; c<save_t.size();c++){
            file<<save_t.at(c)<<"\t";
            for(int c1=0; c1<save_c.at(c).size()-1;c1++) file<<save_c.at(c).at(c1)<<"\t";
            file<<save_c.at(c).at(save_c.at(c).size()-1)<<endl;
        }

        file.close();
    }
    else cout << "Unable to open file './output_curs.dat'"<<endl;
}

void Myocyte::saveExtras(){
    // string curs_path = outputFilePath+"/output_extras.dat";
    // ofstream file (curs_path);
    // if(file.is_open()){
    //     for(int c=0; c<save_t.size();c++){
    //         file<<save_t.at(c)<<"\t";
    //         for(int c1=0; c1<save_e.at(c).size()-1;c1++) file<<save_e.at(c).at(c1)<<"\t";
    //         file<<save_e.at(c).at(save_e.at(c).size()-1)<<endl;
    //     }

    //     file.close();
    // }
    // else cout << "Unable to open file './output_extras.dat'"<<endl;
}

Myocyte::~Myocyte(){
    //for(int c=0; c<tUnits; c++) delete this->ca_units.at(c);
}