strat1 = strat2 = 4 #Initial strategy values: 0.5 for close, 4 for medium, 10 for far
time = 2000 #Time to run simulation
sa=2 #Species niche width that controls adaptive radiation
extinct_fast = 0 #To track extinction times
extinct_slow = 0

IC = [pop1,pop2,strat1,strat2]

#Code to simulate ODEs
def evoLV(X, t):

    gamma = 0

    x1 = X[0]
    x2 = X[1]
    u1 = X[2]
    u2 = X[3]

    K1 = KM * math.exp(-((u1 - gamma) ** 2) / (2 * sk))+1
    K2 = KM * math.exp(-((u2 - gamma) ** 2) / (2 * sk))+1

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

print ('Equilibrium x1: %f' %pop[time][0])
print ('Equilibrium u1: %f' %pop[time][2])

print ('Equilibrium x2: %f' %pop[time][1])
print ('Equilibrium u2: %f' %pop[time][3])

#Obtaining extinction times
for i in range(len(pop[:,1])):
    if (pop[:,1][i])<1:
        print(i)
        print('Fast')
        extinct_fast = i
        break

for j in range(len(pop[:,0])):
    if (pop[:,0][j])<1:
        print(j)
        print('Slow')
        extinct_slow = j
        break
        
if strat1 == 0.5:
    txt = 'Close'
elif strat1 == 4:
    txt = 'Medium'
elif strat1 == 10:
    txt = 'Far'

#Plotting ODE simulation
plt.figure()
plt.subplot(211)
plt.title('Adaptive Radiation: ' + txt)
plt.plot(pop[:,0],label='Slow')
plt.plot(pop[:,1],label='Fast')
if extinct_fast>1:
    #plt.axvline(x=extinct_fast,color='k',ls='--')
    plt.plot(extinct_fast,pop[:,1][extinct_fast],'kD')
elif extinct_slow>1:
    #plt.axvline(x=extinct_slow,color='k',ls='--')
    plt.plot(extinct_slow,pop[:,1][extinct_slow],'kD')
plt.legend()
plt.grid(True)
plt.ylabel('Pop Size, x')
plt.subplot(212)
plt.plot(pop[:,2],label='k = ' + str(k[0]))
plt.plot(pop[:,3],label='k = ' + str(k[1]))
plt.grid(True)
plt.ylabel('Indv Strategy, v')
plt.show()


#Plottng adaptive landscape for a snapshot in time
time_G = 2000 
n=5

fast = []
slow = []

def fastG(v):
    gamma = 0
    x1 = pop[time_G][0]
    x2 = pop[time_G][1]
    u1 = pop[time_G][2]
    u2 = pop[time_G][3]
    a1 = p+(1-p)*(math.exp(-(v-u1)**2/(2*sa)))
    a2 = p+(1-p)*(math.exp(-(v-u2)**2/(2*sa)))
    K2 = KM * math.exp(-((v - gamma) ** 2) / (2 * sk))#+5
    Gfunc2 =r[1]/K2 * (K2 - a1*x1 - a2*x2) - d[1]*k[1]
    return Gfunc2

def slowG(v):
    gamma = 0
    x1 = pop[time_G][0]
    x2 = pop[time_G][1]
    u1 = pop[time_G][2]
    u2 = pop[time_G][3]
    a1 = p+(1-p)*(math.exp(-(v-u1)**2/(2*sa)))
    a2 = p+(1-p)*(math.exp(-(v-u2)**2/(2*sa)))
    K1 = KM * math.exp(-((v - gamma) ** 2) / (2 * sk))#+5
    Gfunc1 = r[0]/K1 * (K1 - a2*x2 - a1*x1) - d[0]*k[0]
    return Gfunc1

for i in np.arange(-n,n,.1):
    fast.append(fastG(i))
    slow.append(slowG(i))

plt.plot(np.arange(-n,n,.1),slow,label='Slow',lw=3)
plt.plot(np.arange(-n,n,.1),fast,label='Fast',lw=3)
plt.plot(pop[time_G][2],slowG(pop[time_G][2]),marker='o',color='c',lw=20,markersize=10)
plt.plot(pop[time_G][3],fastG(pop[time_G][3]),marker='o',color='greenyellow',lw=20,markersize=10)
plt.title('Adaptive Radiation: Time '+str(time_G))
plt.xlabel('Evolutionary Strategy: v')
plt.legend()
#plt.xlim(-1,.1)
plt.ylabel('Fitness: G')
plt.show()
