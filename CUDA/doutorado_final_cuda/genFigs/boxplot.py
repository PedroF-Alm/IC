import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 12})

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
    
    plt.figure()
    plt.boxplot(ints, notch=True, sym='+',showfliers=False,patch_artist=False,showmeans=False)
    plt.ylabel('$\sum^{1s} I_{rel}$')
    plt.xticks([1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19],['1','2','4','8','16','32','64','128','256','512','1024','2048','4096','8192','16384','32768','65536','131072','DET'],rotation=45)
    plt.xlabel('$N_{LCC}$ ($N_{RyR}=5N_{LCC}$)')
    plt.ylim(-0.2,3.7)
    plt.grid(ls=':')
    plt.title(_label_[s]+r' | $I_{rel}\times N_{LCC}$')
    plt.savefig('../exit/'+_scena_[s]+'/irelXnlcc_nO.pdf',bbox_inches='tight')



