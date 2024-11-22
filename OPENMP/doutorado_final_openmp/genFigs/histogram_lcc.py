import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 12})
import numpy as np
import csv
import collections

def getOLCCCurves(scena,n):
    _path_ = '../exit/'+scena+'/'+repr(n)+'_lcc/ca_units/output_extras_u0.dat'
    _olcc_ = 1
    
    var = []
    with open(_path_, 'r') as file_v:
        lines_v = csv.reader (file_v,delimiter='	')
        for row_v in lines_v:
            if(float(row_v[0])>=1000000.):
                var.append(row_v[_olcc_])
    file_v.close()    
    return var

_n_ = 1100
_nlcc_ = [1,2,4,8,16,32,64,128,256,512]#,1024,2048,4096,8192,16384,32768,65536,131072,0]
_scena_ = ['STO_STO','STO_DET','DET_STO']
_label_ = ['$LCC_{STO}$ & $RyR_{STO}$','$LCC_{STO}$ & $RyR_{DET}$','$LCC_{DET}$ & $RyR_{STO}$']

for s in range(len(_scena_)):    
    plt.figure()
    for n in range(len(_nlcc_)):
        data = getOLCCCurves(_scena_[s], _nlcc_[n])

        frequency = np.array(sorted(collections.Counter(data).items()))
        x = [float(x_i) for x_i in frequency[:,0]]
        y = [float(y_i) for y_i in frequency[:,1]]
        
        plt.plot(x, y,label='$N_{LCC}=$'+repr(_nlcc_[n]))
        plt.fill_between(x, y)
       
    plt.xlabel(r'$O_{LCC}$ Opening Fraction (%)')
    plt.ylabel(r'Ocurrence')
    # plt.ylim(0,100)
    plt.yscale('log')
    plt.xlim(-0.05,1.05)
    plt.title(_label_[s])
    plt.legend(loc='right',bbox_to_anchor=(1.35,0.5))
    plt.grid(ls=':')
    plt.savefig('../exit/'+_scena_[s]+'/hist_lcc.pdf',bbox_inches='tight')



