import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 12})
import numpy as np
import csv

def getSplittedCurves(path,col,n):
    curves = []
    var = []
    with open(path, 'r') as file_v:
        lines_v = csv.reader (file_v,delimiter='	')
        for row_v in lines_v:
            data = []
            for i in range (len(row_v)):
                if(i == 0):
                    data.append(float(row_v[i])/1000.)
                else:
                    data.append(float(row_v[i]))
            var.append(data)
    file_v.close()
    var = np.array(var)
    
    for i in range(n):
        curves.append(var[i*1000:(i+1)*1000,col]) #to append the curves
    
    curves = np.array(curves)
    return var[:1000,0],curves
     
#############################################  
_n_ = 1100
_nlcc_ = [1,2,4,8,16,32,64,128,256,512,0]#,1024,2048,4096,8192,16384,32768,65536,131072,0]
_ap_ = 1
_olcc_ = 1
_oryr_ = 2
_irel_ = 3
_cai_ = 3
_scena_ = 'DET_STO'
_label_ = ['$LCC_{DET}$ & $RyR_{STO}$']
_ylim_ = [0,2.25]

for n in _nlcc_:
    plt.figure()
    _varPath_ = '../exit/'+_scena_+'/'+repr(n)+'_lcc//output_vars.dat'
    _curPath_ = '../exit/'+_scena_+'/'+repr(n)+'_lcc/ca_units/output_curs_u0.dat'
    _extPath_ = '../exit/'+_scena_+'/'+repr(n)+'_lcc/ca_units/output_extras_u0.dat'

    path = _extPath_
    col = _cai_
    
    print('Reading simulation: '+path)
    t,curves = getSplittedCurves(path, col, _n_)
    curves = curves*1000.
    for i in range(len(curves)-100):
        plt.plot(t,curves[i+100])
    plt.text(0.8, 0.9*(_ylim_[1]-_ylim_[0])+_ylim_[0], '$N_{LCC}$='+repr(n),bbox=dict(boxstyle="square",ec=(0,0,0),fc=(1,1,1)))
    plt.title(_label_[0])
    plt.grid(ls=':')
    plt.xlim(0,1.)
    plt.ylim(_ylim_)
    plt.ylabel('$[Ca]_i$')
    plt.xlabel('Time (s)')
    plt.savefig('../exit/'+_scena_+'/cai_'+repr(n)+'.pdf',bbox_inches='tight')
    
