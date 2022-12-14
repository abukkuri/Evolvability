import numpy as np
from scipy.integrate import *
import matplotlib.pyplot as plt
import math
from math import exp
import random
import copy

KM = 100 #Carrying capacity
d = 0.05 #Background death rate
r = 0.25 #Intrinsic growth rate
sk = 12.5 #Range of resources
sa = 2 #Species niche width
gamma = 0 #Strategy that maximizes carrying capacity


IC_slow = [] #Initial conditions: slow evolving
IC_fast = [] #Initial conditions: fast evolving
num_slow = 1 #Number of slow evolving morphs
num_fast = 5 #Number of fast evolving morphs
total_num = num_slow+num_fast #Total number of morphs

global schange,fchange1,fchange2,fchange3
schange = False
fchange = False
fchange1 = False
fchange2 = False
fchange3 = False
fchange4 = False
fchange5 = False
fchange6 = False
tryf=[0,0,0,0,0]
trys=[0,0]
p=0.1

fast_spec=[]
slow_spec=[]

time = 20000 #Time for simulation

def gen_slow(n):
    for i in range(n):
        if i==0:
            var1 = 10
            var2 = 0.5
            IC_slow.append(var1)
            IC_slow.append(var2)
        else:
            var1 = 0
            var2 = 0
            IC_slow.append(var1)
            IC_slow.append(var2)
        
def gen_fast(n):
    for i in range(n):
        if i==0:
            var1 = 10
            var2 = 0.5
            IC_fast.append(var1)
            IC_fast.append(var2)
        elif i==1:
            var1 = 0
            var2 = -1.7574494351312362+random.uniform(-.02,.02)
            IC_fast.append(var1)
            IC_fast.append(var2)
        elif i==2:
            var1 = 0
            var2 = -3.1715445689236925+random.uniform(-.02,.02)
            IC_fast.append(var1)
            IC_fast.append(var2)
        elif i==3:
            var1 = 0
            var2 = -1.6769027282461206+random.uniform(-.02,.02)
            IC_fast.append(var1)
            IC_fast.append(var2)
        elif i==4:
            var1 = 0
            var2 = -2.699780467569802+random.uniform(-.02,.02)
            IC_fast.append(var1)
            IC_fast.append(var2)           
            
            
gen_slow(num_slow)
gen_fast(num_fast)

IC = IC_slow+IC_fast

def carry(v):
    return KM * exp(-((v - gamma) ** 2) / (2 * sk))


def evoLV(X, t):

    global fchange,fchange1,fchange2,fchange3,fchange4,fchange5,fchange6, schange

    slows = X[:num_slow*2]
    fasts = X[num_slow*2:]
    
    slow_strats = slows[1::2]
    slow_pops = slows[::2]
    
    fast_strats = fasts[1::2]
    fast_pops = fasts[::2]
    
    temp_slow = []
    temp_fast = []
    
    def comp(n,ss,fs):
        a=0
        i=0
        j=0
        while i < num_slow:
            a += (p+(1-p)*(exp(-(n - ss[i]) ** 2 / (2 * sa))))*slow_pops[i]
            i+=1
        while j<num_fast:
            a += (p+(1-p)*(exp(-(n - fs[j]) ** 2 / (2 * sa))))*fast_pops[j]
            j+=1
        return a
    
    def slow_update():
        i=0
        while i <num_slow:
            G = r/carry(slow_strats[i]) * (carry(slow_strats[i]) - comp(slow_strats[i],slow_strats,fast_strats))-d*0.2
            a = slow_pops[i]*G+trys[i]
            perts = copy.deepcopy(slow_strats)
            perts[i] = slow_strats[i]+.0001 
            Gmod = r/carry(perts[i]) * (carry(perts[i]) - comp(perts[i],perts,fast_strats))-d*0.2
            if slow_pops[i]<.01:
                G=0
                b=0
            else:
                b = 0.2*(Gmod-G)/.0001 
            temp_slow.append(a)
            temp_slow.append(b)
            i+=1
            
    def fast_update():
        i=0
        while i<num_fast:
            G = r/carry(fast_strats[i]) * (carry(fast_strats[i]) - comp(fast_strats[i],slow_strats,fast_strats))-d*0.5
            a = fast_pops[i]*G+tryf[i]
            pertf = copy.deepcopy(fast_strats)
            pertf[i] = fast_strats[i]+.0001 
            Gmod = r/carry(pertf[i]) * (carry(pertf[i]) - comp(pertf[i],slow_strats,pertf))-d*0.5
            if fast_pops[i]<.01:
                G=0
                b=0
            else:
                b = 0.5*(Gmod-G)/.0001 
            temp_fast.append(a)
            temp_fast.append(b)
            i+=1


    #SPECIATION 1

    if abs(fast_strats[0]+1.7574494351312362)<0.02 and t>200 and fchange==False:
        fast_spec.append(t)
        tryf[1]=1
        fchange=True
        
    # if abs(slow_strats[0]-1.3729726690960438)<0.02 and t>200 and schange==False:
    #     print(t,'slow')
    #     trys[1]=1
    #     schange=True
        
    #SPECIATION 2
        
    # if abs(slow_strats[0]-2.5022588316197785)<0.02 and t>1500 and schange==False:
    #     print(t,'slow')
    #     trys[1]=1
    #     schange=True
    if abs(fast_strats[0]+3.1715445689236925)<0.02 and t>1500 and fchange1==False:
        fast_spec.append(t)
        tryf[2]=1
        fchange1=True
    # if abs(fast_strats[1]+0.34385881571564914)<0.02 and t>1500 and fchange2==False:
    #     print(t,'fast1')
    #     fast_spec.append(t)
    #     tryf[2]=1
    #     fchange2=True
    
    #SPECIATION 3
    # if abs(slow_strats[0]-3.5234630444501165)<0.02 and t>4000 and schange==False:
    #     print(t,'slow')
    #     trys[1]=1
    #     schange=True
    # if abs(fast_strats[0]+4.382023230699852)<0.02 and t>4000 and fchange2==False:
    #     print(t,'fast0')
    #     fast_spec.append(t)
    #     tryf[3]=1
    #     fchange2=True
    # if abs(fast_strats[1]-0.8650529169032841)<0.02 and t>4000 and fchange3==False:
    #     print(t,'fast1')
    #     fast_spec.append(t)
    #     tryf[3]=1
    #     fchange3=True
    if abs(fast_strats[2]+1.6769027282461206)<0.02 and t>4000 and fchange2==False:
       fast_spec.append(t)
       tryf[3]=1
       fchange2=True
       
    #SPECIATION 4
    # if abs(slow_strats[0]-4.527229354904877)<0.02 and t>8000 and schange==False:
    #     print(t,'slow')
    #     trys[1]=1
    #     schange=True
    # if abs(fast_strats[0]+5.3667347308562885)<0.02 and t>8000 and fchange3==False:
    #     print(t,'fast0')
    #     fast_spec.append(t)
    #     tryf[4]=1
    #     fchange3=True
    # if abs(fast_strats[1]-1.9841155701887159)<0.02 and t>8000 and fchange4==False:
    #     print(t,'fast1')
    #     fast_spec.append(t)
    #     tryf[4]=1
    #     fchange4=True
    if abs(fast_strats[2]+2.699780467569802)<0.02 and t>8000 and fchange3==False:
        fast_spec.append(t)
        tryf[4]=1
        fchange3=True
    # if abs(fast_strats[3]+0.33246758469704357)<0.02 and t>8000 and fchange6==False:
    #     print(t,'fast3')
    #     fast_spec.append(t)
    #     tryf[4]=1
    #     fchange6=True
        
    slow_update()
    fast_update()
    
    temp = temp_slow+temp_fast
    
    final=np.array(temp)

    return final

intxv = np.array(IC)
time_sp = np.linspace(0,time,time*10-1)
pop = odeint(evoLV, intxv, time_sp,hmax=1)
    
plt.figure()
plt.title('Adaptive Radiation')
j=0
while j<total_num*2:
    if j%2==1:
        if j<num_slow*2:
            print(pop[:,j][-1])
            plt.plot(time_sp,pop[:,j],c='#1f77b4',label='Slow')
        else:
            print(pop[:,j][-1])
            if j==5:
                plt.plot(time_sp[int(fast_spec[0]*10):],pop[int(fast_spec[0]*10):,j],color='#ff7f0e',label='Fast')
            elif j==7:
                plt.plot(time_sp[int(fast_spec[1]*10):],pop[int(fast_spec[1]*10):,j],color='#ff7f0e')
            elif j==9:
                plt.plot(time_sp[int(fast_spec[2]*10):],pop[int(fast_spec[2]*10):,j],color='#ff7f0e')
            elif j==11:
                plt.plot(time_sp[int(fast_spec[3]*10):],pop[int(fast_spec[3]*10):,j],color='#ff7f0e')
            else:
                plt.plot(time_sp,pop[:,j],color='#ff7f0e')
    j+=1
    
plt.grid(True)
plt.legend()
plt.ylabel('Indv Strategy, v')
plt.xlabel('Time')
plt.show()
