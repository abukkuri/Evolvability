import numpy as np
from scipy.integrate import *
import matplotlib.pyplot as plt
import math
from math import exp

#For 2 species

# Model Parameters

pop1 = 10
pop2 = 10
strat1 = strat2 = .5 #.5 is close, 4 is medium, 10 is far

time = 1000
KM = 100
d = [0.05,0.05]
r = [0.25,0.25]
k = [.2,.5]

'''g = np.random.random_sample(time+1)*4-2
g1 = []
for i in range(len(g)):
    if (i%5==0):
        g1.append(g[i])
    else:
        g1.append(g[i-1])'''

IC = [pop1,pop2,strat1,strat2]

#Test different values for this
sk = 12.5
sa = 2 #100 normal, 2 adaptive radiation
B= 0

def evoLV(X, t):

    gamma =  0 #g1[int(t)]

#Evolutionary Rescue
    #AT: Collapse
    '''if t>4800:
        gamma=-24
    elif t>3900:
        gamma = -20
    elif t>3100:
        gamma = -16
    elif t>2400:
        gamma = -12
    elif t>1500: #less: 1800, right: 1500, more: 1200
        gamma = -8
    elif t>800:
        gamma = -4
    else:
        gamma = 0'''

    #AT: IDH
    '''if t>4500:
        gamma = 0
    elif t>3900:
        gamma = -4.5
    elif t>3200:
        gamma = 0
    elif t>2500:
        gamma = -4.5
    elif t>1700:
        gamma = 0
    elif t>800:
        gamma = -4.5
    else:
        gamma = 0'''

    #AT - Collapse (No Aggressiveness)
    '''if t>5000:
        gamma=-24
    elif t>4200:
        gamma = -20
    elif t>3400:
        gamma = -16
    elif t>2500:
        gamma = -12
    elif t>1700:
        gamma = -8
    elif t>900:
        gamma = -4
    else:
        gamma = 0'''

    #AT - IDH No Aggressiveness
    '''if t>4600:
        gamma = 0
    elif t>3950:
        gamma = -4.5
    elif t>3250:
        gamma = 0
    elif t>2550:
        gamma = -4.5
    elif t>1800:
        gamma = 0
    elif t>900:
        gamma = -4.5
    else:
        gamma = 0'''

    #Abbreviated Collapse
    '''if t>1200:
        gamma = -8
    elif t>800:
        gamma = -4
    else:
        gamma = 0'''

    #Abbreviated IDH
    '''if t>1200:
        gamma=0
    elif t>800:
        gamma = -4.5
    else:
        gamma = 0'''

    '''if t>1600:
        gamma=0
    elif t>900:
        gamma = -4.25
    else:
        gamma = 0'''

    #gamma = np.sin(t/50) 

    x1 = X[0]
    x2 = X[1]
    u1 = X[2]
    u2 = X[3]

    #if x1<1:
     #   x1=0
    #if x2<1:
     #   x2=0

    K1 = KM * math.exp(-((u1 - gamma) ** 2) / (2 * sk))
    K2 = KM * math.exp(-((u2 - gamma) ** 2) / (2 * sk))

    a2 = 1 + math.exp(-(u1 - u2 + B) ** 2 / (2 * sa)) - math.exp(-(B ** 2 / (2 * sa)))
    a1 = 1 + math.exp(-(u2 - u1 + B) ** 2 / (2 * sa)) - math.exp(-(B ** 2 / (2 * sa)))

    dx1dt = x1 * (r[0]/K1 * (K1 - a2*x2 - x1) - d[0]*k[0])
    dx2dt = x2 * (r[1]/K2 * (K2 - a1*x1 - x2) - d[1]*k[1])

    dG1dv = r[0]*(-KM*(-2*gamma + 2*u1)*exp(-(-gamma + u1)**2/(2*sk))/(2*sk) + x2*(2*B + 2*u1 - 2*u2)*exp(-(B + u1 - u2)**2/(2*sa))/(2*sa))*exp((-gamma + u1)**2/(2*sk))/KM + r[0]*(-2*gamma + 2*u1)*(KM*exp(-(-gamma + u1)**2/(2*sk)) - x1 - x2*(1 + exp(-(B + u1 - u2)**2/(2*sa)) - exp(-B**2/(2*sa))))*exp((-gamma + u1)**2/(2*sk))/(2*KM*sk)
    dG2dv = r[1]*(-KM*(-2*gamma + 2*u2)*exp(-(-gamma + u2)**2/(2*sk))/(2*sk) + x1*(2*B - 2*u1 + 2*u2)*exp(-(B - u1 + u2)**2/(2*sa))/(2*sa))*exp((-gamma + u2)**2/(2*sk))/KM + r[1]*(-2*gamma + 2*u2)*(KM*exp(-(-gamma + u2)**2/(2*sk)) - x1*(1 + exp(-(B - u1 + u2)**2/(2*sa)) - exp(-B**2/(2*sa))) - x2)*exp((-gamma + u2)**2/(2*sk))/(2*KM*sk)

    dv1dt = k[0] * dG1dv
    dv2dt = k[1] * dG2dv
    dxvdt = np.array([dx1dt, dx2dt, dv1dt, dv2dt])
    return dxvdt

intxv = np.array(IC)
pop = odeint(evoLV, intxv, range(time+1))

#Plotting Adaptive Landscape
time_G = 50

fast = []
slow = []

def fastG(u2):
    gamma = 0
    x1 = pop[time_G][0]
    x2 = pop[time_G][1]
    u1 = pop[time_G][2]
    a1 = 1 + math.exp(-(u2 - u1 + B) ** 2 / (2 * sa)) - math.exp(-(B ** 2 / (2 * sa)))
    K2 = KM * math.exp(-((u2 - gamma) ** 2) / (2 * sk))
    Gfunc2 =r[1]/K2 * (K2 - a1*x1 - x2) - d[1]*k[1]
    return Gfunc2

def slowG(u1):
    gamma = 0
    x1 = pop[time_G][0]
    x2 = pop[time_G][1]
    u2 = pop[time_G][3]
    a2 = 1 + math.exp(-(u1 - u2 + B) ** 2 / (2 * sa)) - math.exp(-(B ** 2 / (2 * sa)))
    K1 = KM * math.exp(-((u1 - gamma) ** 2) / (2 * sk))
    Gfunc1 = r[0]/K1 * (K1 - a2*x2 - x1) - d[0]*k[0]
    return Gfunc1

for i in np.arange(-5,5,.1):
    fast.append(fastG(i))
    slow.append(slowG(i))

plt.plot(np.arange(-5,5,.1),slow,label='Slow')
plt.plot(np.arange(-5,5,.1),fast,label='Fast')
plt.plot(pop[time_G][2],slowG(pop[time_G][2]),marker='o',color='c',lw=20,markersize=5)
plt.plot(pop[time_G][3],fastG(pop[time_G][3]),marker='o',color='greenyellow',lw=20,markersize=5)
plt.title('Adaptive Radiation (No Aggressiveness): Time '+str(time_G))
plt.xlabel('Evolutionary Strategy: v')
plt.legend()
#plt.xlim(-1,.1)
plt.ylabel('Fitness: G')
plt.show()

print ('Equilibrium x1: %f' %pop[time][0])
print ('Equilibrium u1: %f' %pop[time][2])

print ('Equilibrium x2: %f' %pop[time][1])
print ('Equilibrium u2: %f' %pop[time][3])

#Obtaining extinction times
for i in range(len(pop[:,1])):
    if (pop[:,1][i])<1:
        print(i)
        print('Fast')
        break

for j in range(len(pop[:,0])):
    if (pop[:,0][j])<1:
        print(j)
        print('Slow')
        break

plt.figure()
plt.subplot(211)
plt.title('Adaptive Radiation: Close')
plt.plot(pop[:,0],label='Slow')
plt.plot(pop[:,1],label='Fast')
plt.legend()
plt.grid(True)
plt.ylabel('Pop Size, x')
plt.subplot(212)
plt.plot(pop[:,2],label='k = ' + str(k[0]))
plt.plot(pop[:,3],label='k = ' + str(k[1]))
plt.grid(True)
plt.ylabel('Indv Strategy, v')
plt.show()

'''
time_G = 50

fast = []
slow = []

def fastG(v):
    gamma = 0
    x1 = pop[time_G][0]
    x2 = pop[time_G][1]
    u1 = pop[time_G][2]
    u2 = pop[time_G][3]
    a1 = 1 + math.exp(-(v - u1 + B) ** 2 / (2 * sa)) - math.exp(-(B ** 2 / (2 * sa)))
    a2 = 1 + math.exp(-(v - u2 + B) ** 2 / (2 * sa)) - math.exp(-(B ** 2 / (2 * sa)))
    K2 = KM * math.exp(-((v - gamma) ** 2) / (2 * sk))
    Gfunc2 =r[1]/K2 * (K2 - a1*x1 - a2*x2) - d[1]*k[1]
    return Gfunc2

def slowG(v):
    gamma = 0
    x1 = pop[time_G][0]
    x2 = pop[time_G][1]
    u1 = pop[time_G][2]
    u2 = pop[time_G][3]
    a1 = 1 + math.exp(-(v - u1 + B) ** 2 / (2 * sa)) - math.exp(-(B ** 2 / (2 * sa)))
    a2 = 1 + math.exp(-(v - u2 + B) ** 2 / (2 * sa)) - math.exp(-(B ** 2 / (2 * sa)))
    K1 = KM * math.exp(-((v - gamma) ** 2) / (2 * sk))
    Gfunc1 = r[0]/K1 * (K1 - a2*x2 - a1*x1) - d[0]*k[0]
    return Gfunc1

for i in np.arange(-5,5,.1):
    fast.append(fastG(i))
    slow.append(slowG(i))

plt.plot(np.arange(-5,5,.1),slow,label='Slow')
plt.plot(np.arange(-5,5,.1),fast,label='Fast')
plt.plot(pop[time_G][2],slowG(pop[time_G][2]),marker='o',color='c',lw=20,markersize=5)
plt.plot(pop[time_G][3],fastG(pop[time_G][3]),marker='o',color='greenyellow',lw=20,markersize=5)
plt.title('Adaptive Radiation (No Aggressiveness): Time '+str(time_G))
plt.xlabel('Evolutionary Strategy: v')
plt.legend()
#plt.xlim(-1,.1)
plt.ylabel('Fitness: G')
plt.show()
'''
