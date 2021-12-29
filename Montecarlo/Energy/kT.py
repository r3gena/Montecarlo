import numpy as np

T=1.353113
for i in range(10):
    energy_file="energy_"+str(i+1)+".txt"
    data=np.loadtxt(energy_file).T
    mdphi=np.sum(data[1])/np.size(data[1])
    md2phi=np.sum(data[2])/np.size(data[2])
    mdphi2=np.sum((data[1]*data[1]))/np.size(data[1])
    inv_kT=(500*T/1000)+1000*md2phi-1000*(mdphi2-(mdphi*mdphi))/T
    kT=1/inv_kT
    print("Proceso "+str(i+1)+": ",kT)
