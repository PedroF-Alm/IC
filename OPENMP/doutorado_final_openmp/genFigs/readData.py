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

def getSplittedCurvesIntegral(path,col,n):
    t,curves = getSplittedCurves(path, col, n)
    ints = []
    for c in curves[100:]:
        ints.append(sum(c))
    ints = np.array(ints)
    return ints
        
#############################################  
_n_ = 1100
_nlcc_ = [1,2,4,8,16,32,64,128,256,512,1024,2048,4096,8192,16384,32768,65536,131072,0]
_irel_ = 3
_ical_ = 5
_olcc_ = 1
_oryr_ = 2

ints = []

for n in _nlcc_:
    _curPath_ = '../exit/STO_STO/'+repr(n)+'_lcc/ca_units/output_curs_u0.dat'
    _extPath_ = '../exit/STO_STO/'+repr(n)+'_lcc/ca_units/output_extras_u0.dat'

    path = _curPath_
    col = _irel_
    print('Reading simulation: '+path)
    ints.append(getSplittedCurvesIntegral(path,col,_n_))

plt.figure()
plt.boxplot(ints, notch=True, sym='+',showfliers=False)
plt.ylabel('$\sum^{1s} I_{rel}$')
plt.xticks([1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19],['1','2','4','8','16','32','64','128','256','512','1024','2048','4096','8192','16384','32768','65536','131072','DET'],rotation=45)
plt.xlabel('$N_{LCC}$ ($N_{RyR}=5N_{LCC}$)')
plt.grid(ls=':')
plt.title(r'$I_{rel}\times N_{LCC}$ (n=1000)')
plt.savefig('irelXnlcc.pdf',bbox_inches='tight')



