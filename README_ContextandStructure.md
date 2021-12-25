# Context and Structure

## Physical Framework
The code presented here is part of my thesis work. The aim is to create a simple one-dimensional model that simulate the beam diffusion inside a circular accelerator and to compare the result with the solution of a Fokker-Planck equation. More in detail we study the use of a relatively new tool to implement in HL-LHC, namely an Hollow Electron Lens (HEL), in order to boost the performance of the collimation system through an active control of halo particles’ diffusion speed and tail population.
The lens is thus an hollow beam of electrons that acts only on the halo of the protons beam at a specific point of the machine, while leaving the core unchanged. This translate in to a kick in the map of the protons that can be modulated and made stochastic in order to create the necessary noise to produce the diffusion. The final goal would be to find the noise that can allow the fastest diffusion possible.

In this code we generate initial disrtibutions of protons and we choose various simulation's parameters as the HEL inner and outer radius (given in RMS beam size), or the number of initial protons and iterations of the map. All the other parameters are given by the LHC working conditions. It is possible to change the kind of noise to have different results; in the module *core.py*, inside the *make_correlated_noise* function, few noise examples are commented. There can also be a correlation in the noise realizations which depends on the *gamma* parameter (dafault is 0). 
It is possible to write a script that ask to run the full simulation and print results, as the one here, or to run the simulation by steps and see the gradual evolution of the system. 
After the simulation the results outputs are plotted and saved, and the information about the halo loss percentage is printed.

## Organization of the code
In the **Project folder** one can find the codes and output files:

 --------------------------------------------------------------- **Code Files** ---------------------------------------------------------------
 
 1. *script_map.py*: Script in which initial distributions are generated and parameteres are set. Here an instance of the simulation is generated and the iterations are called. Then, all the results are saved and plotted;

 2. *_initmap_.py*: Module containing the classes, with various useful methods to monitor the simulation, including the one that calculate the current of particles loss.

 3. *core.py*: Module containing all the actual mathematical calculations of the map iterations.

 --------------------------------------------------------------- **Output Files** --------------------------------------------------------------
 
 9. *data_I_x_p_th.txt*: Output file contianing data of the final action, position, momentum and angle of the surviving halo particles;
 10. *Action_distribution.pdf*: Pdf file containing the plot of the final action distribution;
 11. *Theta_distribution.pdf*: Pdf file containing the plot of the final angle distribution;
 12. *Current.pdf*: Pdf file containing the plot of the particle loss current;
 13. *Phase_space.pdf*: Pdf file containing the plot of the final halo phase space.
