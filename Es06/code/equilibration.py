import matplotlib
import matplotlib.pyplot as plt
import numpy as np

def loadFile(name):
    dataFile = open(name, "r")
    lines = dataFile.readlines()
    fileLength = len(lines)
    data = np.empty(fileLength)
    for i in range(fileLength):
        values = lines[i].split()
        data[i] = float(values[1])

    dataFile.close()

    return data, fileLength



path = "gibbs/"

fig, ax = plt.subplots(2, 2, figsize = (12, 8))

points = 100
T = 2.0
beta = 1/T
J = 1.0
Ns = 50
th = np.tanh(J/T)
thN= th**Ns
ch = 1/th


data, fileLength = loadFile(path + "output.ene.0")
x = np.arange(fileLength) + 1

ax[0, 0].plot(x, data)
ax[0, 0].set(xlabel = "#step", ylabel = "U/N")
ax[0, 0].grid(True)

e = -J*( th + ch*thN )/( 1 + thN )
ax[0, 0].axhline(e, color = "#d82786")


data, fileLength = loadFile(path + "output.heat.0")
x = np.arange(fileLength) + 1

ax[0, 1].plot(x, data)
ax[0, 1].set(xlabel = "#step", ylabel = "C")
ax[0, 1].grid(True)

heat=((beta*J)**2)*(((1+thN+(Ns-1)*(th**2)+(Ns-1)*(ch**2)*thN)/(1+thN))-Ns*((th+ch*thN)/(1+thN))**2)
ax[0, 1].axhline(heat, color = "#d82786")


data, fileLength = loadFile(path + "output.mag.0")
x = np.arange(fileLength) + 1

ax[1, 0].plot(x, data)
ax[1, 0].set(xlabel = "#step", ylabel = "M")
ax[1, 0].grid(True)

h=0.02 #external field

l1 = np.exp(beta*J)*np.cosh(beta*h)+np.sqrt(np.exp(2*beta*J)*np.cosh(beta*h)*np.cosh(beta*h)-2*np.sinh(2*beta*J))
l2 = np.exp(beta*J)*np.cosh(beta*h)-np.sqrt(np.exp(2*beta*J)*np.cosh(beta*h)*np.cosh(beta*h)-2*np.sinh(2*beta*J))
Z = l1**Ns + l2**Ns
M = (np.exp(beta*J)*np.sinh(beta*h)*((l1**(Ns-1))*(1+np.exp(beta*J)*np.cosh(beta*h)/np.sqrt(np.exp(2*beta*J)*np.cosh(beta*h)*np.cosh(beta*h)-2*np.sinh(2*beta*J))) 
        + (l2**(Ns-1))*(1-np.exp(beta*J)*np.cosh(beta*h)/np.sqrt(np.exp(2*beta*J)*np.cosh(beta*h)*np.cosh(beta*h)-2*np.sinh(2*beta*J)))))/(Z)
ax[1, 0].axhline(M, color = "#d82786")


data, fileLength = loadFile(path + "output.chi.0")
x = np.arange(fileLength) + 1

ax[1, 1].plot(x, data)
ax[1, 1].set(xlabel = "#step", ylabel = "$\chi$")
ax[1, 1].grid(True)

X = beta*np.exp(2*beta*J)*(1-thN)/(1+thN)
ax[1, 1].axhline(X, color = "#d82786")


fig.tight_layout()
plt.show()