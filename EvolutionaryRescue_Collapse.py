strat1 = strat2 = 4
time = 5000
sa=100
extinct_fast = 0
extinct_slow = 0

IC = [pop1,pop2,strat1,strat2]

#Code to simulate ODE
def evoLV(X, t):
    
    if t>4200:
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
        gamma = 0

    x1 = X[0]
    x2 = X[1]
    u1 = X[2]
    u2 = X[3]
    
    if x1<1:
        x1=0
    if x2<1:
        x2=0

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
        
if dist == 1900:
    txt = 'Rare'
elif dist == 1600:
    txt = 'Just Right'
elif dist == 1400:
    txt = 'Frequent'

#Code to plot ODE
plt.figure()
plt.subplot(211)
plt.title('Cyclical Coexistence: Ecological Collapse')
plt.plot(pop[:,0],label='Slow')
plt.plot(pop[:,1],label='Fast')
if extinct_fast>1:
    plt.plot(extinct_fast,pop[:,1][extinct_fast],'kD')
elif extinct_slow>1:
    plt.plot(extinct_slow,pop[:,0][extinct_slow],'kD')
plt.legend()
plt.grid(True)
plt.ylabel('Pop Size, x')
plt.subplot(212)
plt.plot(pop[:,2],label='k = ' + str(k[0]))
plt.plot(pop[:,3],label='k = ' + str(k[1]))
plt.grid(True)
plt.ylabel('Indv Strategy, v')
plt.show()
