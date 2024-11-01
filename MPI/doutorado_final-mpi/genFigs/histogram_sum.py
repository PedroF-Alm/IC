import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 12})
import numpy as np

_n_ = 1100
_nlcc_ = [1,2,4,8,16,32,64,128,256,512,1024,2048,4096,8192,16384,32768,65536,131072,0]
_scena_ = ['STO_STO','STO_DET','DET_STO']
_label_ = ['$LCC_{STO}$ & $RyR_{STO}$','$LCC_{STO}$ & $RyR_{DET}$','$LCC_{DET}$ & $RyR_{STO}$']

for s in range(len(_scena_)):
    ints = []
    for n in _nlcc_:
        file = '../exit/'+_scena_[s]+'/'+repr(n)+'_lcc/irelSum.dat'
    
        print('Reading simulation: '+file)
        data = []
        file_v = open(file, 'r')
        for l in file_v:
            data.append(float(l))
        file_v.close()
        ints.append(data)
    
    for n in range(len(_nlcc_)):
        plt.figure()
        bins_list = np.arange(0,3.5,0.1)
        plt.hist(ints[n],bins=bins_list,color = 'red', edgecolor = 'black',label='$N_{LCC}$='+repr(_nlcc_[n]))
        plt.xlabel('$\sum^{1s} I_{rel}$')
        plt.ylabel('Ocurrence')
        plt.ylim(0,500)
        plt.xlim(-0.2,3.7)
        plt.grid(ls=':')
        plt.title(_label_[s])
        plt.legend()
        plt.savefig('../exit/'+_scena_[s]+'/hist_'+repr(_nlcc_[n])+'.pdf',bbox_inches='tight')



