import numpy as np
from numba import njit, prange


@njit
def f (x, R_1, R_2):
    """internal calculation of the kick function

    Parameters
    ----------
    x : float
        position
    R_1 : float
        inner lens radius
    R_2 : float
        outer lens radius

    Returns
    -------
    float 
        kick value
    """

    if abs(x) < R_1:
        return 0.0
        
    if abs(x) > R_2:
        return 1.0

    else :
        result= (x**2 - R_1**2)/(R_2**2 - R_1**2)
        return result


def make_correlated_noise(n_elements, gamma=0.0):
    """Make an array of correlated noise
    
    Parameters
    ----------
    n_elements : unsigned int
        number of elements
    gamma : float, optional
        correlation coefficient, by default 0.0
    
    Returns
    -------
    ndarray
        the noise array
    """
    #one can choose any kind of noise, here are listed few examples
    #noise=np.full(n_elements, 0.5)
    #noise=np.random.rand(n_elements)
    #noise = np.random.normal(0.0, 1.0, n_elements)
    noise = np.random.binomial(1, 0.5, n_elements)
    if gamma != 0.0:
        for i in range(1, n_elements):
            noise[i] += gamma * noise[i - 1]
    return noise


@njit
def iterate(x, px, noise, omega_0, omega_1, omega_2, R_1, R_2, TH_MAX, barrier_radius, start):
    """internal iteration method for symplectic map
    
    Parameters
    ----------
    x : float
        x0
    px : float
        px0
    noise : ndarray
        array of noise values
    omega_0 : float
        omega 0
    omega_1 : float
        omega 1
    omega_2 : float
        omega 2
    R_1 : float
        inner lens radius
    R_2 : float
        outer lens radius
    TH_MAX : float
        theta max value
    barrier_radius : float
        barrier position 
    start : unsigned int
        starting iteration value
    
    Returns
    -------
    (float, float, unsigned int)
        (x, px, iterations)
    """    
    for i in range(len(noise)):        
        action = (x * x + px * px) * 0.5
        rot_angle = omega_0 + (omega_1 * action) + (0.5*  omega_2 * action * action) 
        if (x==0) and (px==0):
            return 0.0, 0.0, start
        if (np.sqrt(action*2) >= barrier_radius):
            return 0.0, 0.0, i + start 
        temp1 = x
        temp2 = (px + (noise[i]
         * TH_MAX * f(x, R_1, R_2) * R_2 * 914099.8357243269 * x/ (abs(x)*abs(x))))
        x = np.cos(rot_angle) * temp1 + np.sin(rot_angle) * temp2
        px = -np.sin(rot_angle) * temp1 + np.cos(rot_angle) * temp2
        action = (x * x + px * px) * 0.5
        if (np.sqrt(action*2) >= barrier_radius):
            return 0.0, 0.0, i + start 
        
    return x, px, i + start +1



@njit(parallel=True)
def symplectic_map_common(x, px, step_values, noise_array, omega_0, omega_1, omega_2, R_1, R_2, TH_MAX, barrier_radius, gamma=0.0):
    """computation for common noise symplectic map

    Parameters
    ----------
    x : ndarray
        x initial condition
    px : ndarray
        px initial condition
    step_values : ndarray
        iterations already performed
    noise_array : ndarray
        noise array for the whole group
    omega_0 : float
        omega 0
    omega_1 : float
        omega 1
    omega_2 : float
        omega 2
    R_1 : float
        inner lens radius
    R_2 : float
        outer lens radius
    TH_MAX : float
        theta max value
    barrier_radius : float
        barrier radius
    gamma : float, optional
        correlation coefficient, by default 0.0

    Returns
    -------
    (ndarray, ndarray, ndarray)
        x, px, step_values
    """
    for i in prange(len(x)):
        x[i], px[i], step_values[i] = iterate(x[i], px[i], noise_array, omega_0, omega_1, omega_2, R_1, R_2, TH_MAX, barrier_radius,  step_values[i])
    return x, px, step_values 
