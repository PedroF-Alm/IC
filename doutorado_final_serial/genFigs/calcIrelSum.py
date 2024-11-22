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
_scena_ = 'DET_STO'

ints = []

for n in _nlcc_:
    _curPath_ = '../exit/'+_scena_+'/'+repr(n)+'_lcc/ca_units/output_curs_u0.dat'

    path = _curPath_
    col = _irel_
    print('Reading simulation: '+path)
    ints = getSplittedCurvesIntegral(path,col,_n_)
    outputFile = open('../exit/'+_scena_+'/'+repr(n)+'_lcc/irelSum.dat','w')
    for i in range(len(ints)):
        outputFile.write(repr(ints[i])+'\n')



