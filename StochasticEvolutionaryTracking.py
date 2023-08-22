import numpy as np
from scipy.integrate import *
import matplotlib.pyplot as plt
import math
from math import exp

pop1 = 10
pop2 = 10
strat1 = strat2 = .5

time = 1000
KM = 100
d = [0.05,0.05]
r = [0.25,0.25]
k = [.2,.5] 

IC = [pop1,pop2,strat1,strat2]

sk = 12.5
sa = 100
B = 0

pop_macro = []

#Code for ODE simulation: one run
def run_this():

    g = np.random.random_sample(time + 1) * 4 - 2
    g1 = []
    for i in range(len(g)):
        if (i % 5 == 0):
            g1.append(g[i]) 
        else:
            g1.append(g[i - 1])

    def evoLV(X, t):
        gamma =  g1[int(t)]

        x1 = X[0]
        x2 = X[1]
        u1 = X[2]
        u2 = X[3]
        
        if x1<2:
            x1=0
            dx1dt=0
        if x2<2:
            x2=0
            dx2dt=0

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
    pop_macro.append(pop)

num_sims = 100 #number of simulation runs
for i in range(0,num_sims):
    run_this()

fast_extinct = []
slow_extinct = []

#Obtaining extinction times
for j in range(0,num_sims):
    for i in range(len(pop_macro[j][:,1])):
        if (pop_macro[j][:,1][i])<2:
            fast_extinct.append(i)
            break

for j in range(0,num_sims):
    for i in range(len(pop_macro[j][:,0])):
        if (pop_macro[j][:,0][i])<2:
            slow_extinct.append(i)
            break

line_width = .3

#Plot stochastic ODE simulations
plt.figure()
plt.subplot(211)
plt.title('Stochastic Evolutionary Tracking: Close')
plt.plot(pop_macro[0][:,0],label='Slow',lw=1,c='C1')
plt.plot(pop_macro[0][:,1],label='Fast',lw=1,c='C0')
plt.plot(fast_extinct[0],pop_macro[0][:,1][fast_extinct[0]],'kD')
plt.plot(slow_extinct[0],pop_macro[0][:,0][slow_extinct[0]],'kD')
for i in range(1,num_sims):
    plt.plot(pop_macro[i][:,0],lw=line_width,c='C1')
    plt.plot(pop_macro[i][:,1],lw=line_width,c='C0')
    plt.plot(fast_extinct[i],pop_macro[i][:,1][fast_extinct[i]],'kD')
    plt.plot(slow_extinct[i],pop_macro[i][:,0][slow_extinct[i]],'kD')
plt.legend()
plt.grid(True)
plt.ylabel('Pop Size, x')
plt.subplot(212)
plt.plot(pop_macro[0][:,2],lw=line_width,c='C1')
plt.plot(pop_macro[0][:,3],lw=line_width,c='C0')
for i in range(1,num_sims):
    plt.plot(pop_macro[i][:,2],lw=line_width,c='C1')
    plt.plot(pop_macro[i][:,3],lw=line_width,c='C0')
plt.grid(True)
plt.ylabel('Indv Strategy, v')
plt.show()
