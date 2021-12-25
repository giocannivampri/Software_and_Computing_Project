from init_map import make_correlated_noise, symplectic_map as sm
import numpy as np
import matplotlib.pyplot as plt
np.random.seed(4)


#machine's parameters (these are given and fixed based on the machine's condition)
radius=6.7
th_MAX=  0.3e-6
emittance_star=2.5e-6
beta=280
gamma_rel=(7000/0.938272088)
beta_rel=np.sqrt(1-(1/(gamma_rel**2)))
emittance=emittance_star/(gamma_rel*beta_rel)
omega_0x=0.31 * 2 * np.pi
omega_1x= -1.73e5 * 2*np.pi * 2*emittance_star/(beta_rel*gamma_rel)
omega_2x= -1.87e12 * 2*np.pi * (2*emittance_star/(beta_rel*gamma_rel))**2
mean_x=0.0
sigma_x=1.0
sigma_x2=2.0
mean_p=0.0
sigma_p=1.0
sigma_p2=2.0


#setting simulation's parameters (one can play around with these to generate various results)
r_1=5   #inner radius of the lens expressed in RMS beam size
r_2=2*r_1   #outer radius of the lens, usually 2 times r_1
iterations= 10**5   #iterations of the simulation, i.e. the number of turns in the accelerator
n_particle=10**6   #number of total initial beam protons



#creating initial distributions
Min_I=r_1**2 * 0.5 
Max_I=radius**2 * 0.5
x0=np.random.normal(mean_x, sigma_x, int(n_particle*0.65))
x1=np.random.normal(mean_x, sigma_x2, int(n_particle*0.35))
x0=np.append(x0, x1)
p0=np.random.normal(mean_p, sigma_p, int(n_particle*0.65))
p1=np.random.normal(mean_p, sigma_p2, int(n_particle*0.35))
p0=np.append(p0, p1)
I0=(x0**2 + p0**2)*0.5
p0=(p0)[I0>Min_I]
x0=(x0)[I0>Min_I]
I0=(I0)[I0>Min_I]
p0=(p0)[I0<Max_I]
x0=(x0)[I0<Max_I]
I0=(I0)[I0<Max_I]


#generating instance and iterations of the map
mappa = sm.generate_instance(omega_0x, omega_1x, omega_2x, r_1, r_2, th_MAX, radius, x0, p0)
mappa.common_noise(make_correlated_noise(iterations))


#output
th=mappa.get_th()
y=mappa.get_filtered_action()
x, p, t = mappa.get_filtered_data()
per=(len(x0)-len(x))/len(x0)
print('__________________________________________________\n')
print('Percentage of beam affected by HEL (Halo):', len(x0)/n_particle)
print('Halo loss percentage after', iterations,'iterations:',per)
print('__________________________________________________\n')
binning=3000
tempi_mappa, corrente_mappa = mappa.current_binning(binning)

with open('output/data_I_x_p_th.txt', 'w') as f:
    f.write('')
    f.close()

with open('output/data_I_x_p_th.txt', 'a') as f:
    for i in range(len(y)):
        f.write(str(y[i]))
        f.write(' ')
        f.write(str(x[i]))
        f.write(' ')
        f.write(str(p[i]))
        f.write(' ')
        f.write(str(th[i]))
        f.write('\n')
    f.close()

fig, d=plt.subplots()
d.set_xlabel("I")
d.set_ylabel("\u03C1")
d.hist(I0, 30, density=True, label='Initial distribution')
d.hist(y, 30, density=True, histtype= 'step', label='Final distribution')
plt.legend(loc='best')
plt.savefig('output/Action_distribution.pdf')

fig, d=plt.subplots()
d.set_xlabel("\u03B8")
d.set_ylabel("\u03C1")
d.hist(th, 50, density=True)
plt.savefig('output/Theta_distribution.pdf')

fig, d=plt.subplots()
d.set_xlabel("Turns")
d.set_ylabel("Current (particles/turn)")
d.plot(tempi_mappa, corrente_mappa)
plt.savefig('output/Current.pdf')

fig, d=plt.subplots()
d.set_xlabel("x")
d.set_ylabel("P")
d.plot(x, p, "r.", markersize=1)
plt.savefig('output/Phase_space.pdf')

plt.show()
