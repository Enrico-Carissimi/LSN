import matplotlib
import matplotlib.pyplot as plt
import numpy as np

def loadFile(name):
    dataFile = open(name, "r")
    lines = dataFile.readlines()
    fileLength = len(lines)
    temp = np.empty(fileLength)
    mean = np.empty(fileLength)
    uncertainty = np.empty(fileLength)
    
    for i in range(fileLength):
        values = lines[i].split()
        temp[i], mean[i], uncertainty[i] = float(values[0]), float(values[1]), float(values[2])

    dataFile.close()

    return temp, mean, uncertainty, fileLength



fig, ax = plt.subplots(2, 2, figsize = (12, 8))


path = "gibbs/"

points=100
T = np.linspace(0.4, 2.1, num=points)
beta = 1/T
J = 1.0
Ns = 50
th = np.tanh(J/T)
thN= th**Ns
ch = 1/th


#energy
temp, mean, uncertainty, fileLength = loadFile(path + "output.ene.final")

ax[0, 0].errorbar(temp, mean, yerr = uncertainty, fmt = "o", color = "black", ecolor = "#1f77b4")
ax[0, 0].set(xlabel = "#blocks", ylabel = "U/N")
ax[0, 0].grid(True)

e = -J*( th + ch*thN )/( 1 + thN )
ax[0, 0].plot(T, e, color = "#d82786")


#heat capacity
temp, mean, uncertainty, fileLength = loadFile(path + "output.heat.final")
x = np.arange(fileLength) + 1

ax[0, 1].errorbar(temp, mean, yerr = uncertainty, fmt = "o", color = "black", ecolor = "#1f77b4")
ax[0, 1].set(xlabel = "#blocks", ylabel = "C []")
ax[0, 1].grid(True)

heat=((beta*J)**2)*(((1+thN+(Ns-1)*(th**2)+(Ns-1)*(ch**2)*thN)/(1+thN))-Ns*((th+ch*thN)/(1+thN))**2)
ax[0, 1].plot(T, heat, color = "#d82786")


#magnetization
temp, mean, uncertainty, fileLength = loadFile(path + "output.mag.final")
x = np.arange(fileLength) + 1

ax[1, 0].errorbar(temp, mean, yerr = uncertainty, fmt = "o", color = "black", ecolor = "#1f77b4")
ax[1, 0].set(xlabel = "#blocks", ylabel = "M []")
ax[1, 0].grid(True)

h=0.02 #external field

l1 = np.exp(beta*J)*np.cosh(beta*h)+np.sqrt(np.exp(2*beta*J)*np.cosh(beta*h)*np.cosh(beta*h)-2*np.sinh(2*beta*J))
l2 = np.exp(beta*J)*np.cosh(beta*h)-np.sqrt(np.exp(2*beta*J)*np.cosh(beta*h)*np.cosh(beta*h)-2*np.sinh(2*beta*J))
Z = l1**Ns + l2**Ns
M = (np.exp(beta*J)*np.sinh(beta*h)*((l1**(Ns-1))*(1+np.exp(beta*J)*np.cosh(beta*h)/np.sqrt(np.exp(2*beta*J)*np.cosh(beta*h)*np.cosh(beta*h)-2*np.sinh(2*beta*J))) 
        + (l2**(Ns-1))*(1-np.exp(beta*J)*np.cosh(beta*h)/np.sqrt(np.exp(2*beta*J)*np.cosh(beta*h)*np.cosh(beta*h)-2*np.sinh(2*beta*J)))))/(Z)
ax[1, 0].plot(T, M, color = "#d82786")


#magnetic susceptibility
temp, mean, uncertainty, fileLength = loadFile(path + "output.chi.final")
x = np.arange(fileLength) + 1

ax[1, 1].errorbar(temp, mean, yerr = uncertainty, fmt = ".", color = "black", ecolor = "#1f77b4")
ax[1, 1].set(xlabel = "#blocks", ylabel = "chi []")
ax[1, 1].grid(True)

X = beta*np.exp(2*beta*J)*(1-thN)/(1+thN)
ax[1, 1].plot(T, X, color = "#d82786")


fig.tight_layout()
plt.show()