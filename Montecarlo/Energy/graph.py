import numpy as np
import matplotlib.pyplot as plt

Ep=np.array([])
dphi=np.array([])
d2phi=np.array([])
for i in range(10):
    energy_file="energy_"+str(i+1)+".txt"
    data=np.loadtxt(energy_file).T
    Ep=np.append(Ep,data[0])
    dphi=np.append(dphi,data[1])
    d2phi=np.append(d2phi,data[2])
plt.figure(1)
plt.plot(Ep)
plt.title("Potencial")
plt.ylabel("Ep")
plt.savefig("Ep.png")

plt.figure(2)
plt.plot(dphi)
plt.title("Primera derivada")
plt.ylabel("dphi")
plt.savefig("dphi.png")

plt.figure(3)
plt.plot(d2phi)
plt.title("Segunda derivada")
plt.ylabel("d2phi")
plt.savefig("d2phi.png")
