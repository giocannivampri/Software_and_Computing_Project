# Optimizations
This code have to compute many iterations of the simulation and it does that for many initial conditions. If we consider for example a distribution of 10^6 initial particle and 10^5 iterations of the simulation, we would need around 10^11 total iteration calculations. Of course the first trivial way to lighten the simulation is to simulate only the halo, which is the part we are interested in, given the fact that in the core the HEL field is null. This translate in a condition on the action variable of the particles; so we simulate only the particle for which the action is larger than $r_1^2/2$, where r_1 is the HEL inner radius.
In this way we are left with only around $10^4$ particles to iterate 10^5 times. It is still a large amount of computation to make.
To speed up the simulation here we make use of Numba. Numba is a just-in-time compiler, really useful if one have code that uses NumPy arrays and functions with loops. When a call is made to a Numba-decorated function, it is compiled to machine code “just-in-time” for execution and the code can run at native machine code speed. 

The runtime of the simulation of $10^6$ total initial particles (around $14*10^3$ halo particles fixing r_1 to 5) for 10^3 iterations is around 198 seconds, while using Numba it is around 1.4 senconds, which gives a **speed up** of around 141.4.

A further optimization is the the option to enable the automatic parallelization of the code using Numba by specifying *Parallel=True* in the decorator. If, thiss time, we perform $2*10^5$ iterations, the regular simulation takes around 70 seconds, while the parallelized one takes around 24 seconds, which gives a further **speed up** of around 2.9.


