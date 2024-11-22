import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 12})
import numpy as np

_x_ = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19]
_nlcc_ = [1,2,4,8,16,32,64,128,256,512,1024,2048,4096,8192,16384,32768,65536,131072,0]
_scena_ = ['STO_STO','STO_DET','DET_STO']
_label_ = ['$LCC_{STO}$ & $RyR_{STO}$','$LCC_{STO}$ & $RyR_{DET}$','$LCC_{DET}$ & $RyR_{STO}$']

for s in range(len(_scena_)):
    avg = []
    std = []
    his = []
    for n in _nlcc_:
        file = '../exit/'+_scena_[s]+'/'+repr(n)+'_lcc/irelSum.dat'
    
        print('Reading simulation: '+file)
        data = []
        file_v = open(file, 'r')
        for l in file_v:
            data.append(float(l))
        file_v.close()
        data = np.sort(data)
        his.append(data[-10:])
        avg.append(np.mean([data]))
        std.append(np.std([data]))
    
    plt.figure()
    plt.errorbar(_x_,avg,std,color='red',marker='o',label='MEAN$\pm$STD')
    for l in range(len(_nlcc_)):
        aux = 10*[_x_[l]]
        plt.scatter(aux,his[l],marker='+',color='black')
    plt.scatter(None,None,marker='+',color='black',label='10 Highest')
    plt.ylabel('$\sum^{1s} I_{rel}$')
    plt.xticks(_x_,['1','2','4','8','16','32','64','128','256','512','1024','2048','4096','8192','16384','32768','65536','131072','DET'],rotation=45)
    plt.xlabel('$N_{LCC}$ ($N_{RyR}=5N_{LCC}$)')
    plt.ylim(-0.2,3.7)
    plt.grid(ls=':')
    plt.legend()
    plt.title(_label_[s]+r' | $I_{rel}\times N_{LCC}$')
    plt.savefig('../exit/'+_scena_[s]+'/irelXnlcc_avgStd.pdf',bbox_inches='tight')



