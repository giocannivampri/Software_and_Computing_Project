# How To Run
## 0) Install Numba
if you haven't already installed it, Numba is available as a conda package for the Anaconda Python distribution. Once you have conda installed, just type:

    $ conda install numba

or you can install Numba using pip:

    $ pip install numba

## 1) Download and unzip

Download the project from the repository of GitHub and unzip the folder to extract the files.

## 2) Run the code

Enter on the bash the position of the downloaded folder:

    cd path/Software_and_Computing_Project-main/project

Then type:

    python script_map.py

All the outputs are saved in the *output* folder


## 3) Play around with the simulation

It is possible to personalize the simulation by modifying the parameters inside *script.py*; one can change the number of initial particles and the number of iterations of the map, and it is also possible to choose the HEL radius. One can also experiment with variuos noises by changing it inside the *make_correlated_noise* function present in *core.py* (a few examples of noise are commented here), or generate a correlated noise by choosing a correlation coefficient (*gamma*) while calling the aforementioned function.
