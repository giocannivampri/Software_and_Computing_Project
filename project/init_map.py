import numpy as np
import core as mio

def make_correlated_noise(n_elements, gamma=0.0):
    return mio.make_correlated_noise(n_elements, gamma)

class symplectic_map(object):
    def __init__(self):
        pass
    
    def reset(self):
        pass

    def common_noise(self):
        pass
    
    def get_data(self):
        """Get data from engine.
        
        Returns
        -------
        tuple(ndarray, ndarray, ndarray)
            tuple with x, p, and number of iterations before loss data
        """
        return self.x, self.p, self.times

    def get_filtered_data(self):
        """Get filtered data from engine.
        
        Returns
        -------
        tuple(ndarray, ndarray, ndarray)
            tuple with x, p, and number of iterations before loss data
        """
        t = self.times
        return (self.x)[t >= self.iterations], (self.p)[t >= self.iterations], t[t >= self.iterations]
   
    def get_th(self):
        """Find angle variable data from engine.
        
        Returns
        -------
        ndarray
            angle array data
        """
        x, p, t = self.get_filtered_data()
        th=[]
        for i in range(len(x)):
            th_1=np.arcsin(x[i]/np.sqrt((x[i] * x[i])+(p[i] * p[i])))
            th_2=np.arccos(p[i]/np.sqrt((x[i] * x[i])+(p[i] * p[i])))
            if np.sin(th_1) > 0 and np.cos(th_2) > 0 :
                th.append(th_1)
            if np.sin(th_1) > 0 and np.cos(th_2) < 0 :
                th.append(th_2)
            if np.sin(th_1) < 0 and np.cos(th_2) > 0 :
                th_4=(th_1 + (np.pi*2))
                th.append(th_4)
            if np.sin(th_1) < 0 and np.cos(th_2) < 0 :
             th_3=(np.pi - th_1)
             th.append(th_3)
            
        Th=np.array(th)
        return Th

    def get_action(self):
        """Get action data from engine
        
        Returns
        -------
        ndarray
            action array data
        """
        return (self.x * self.x + self.p * self.p) * 0.5

    def get_filtered_action(self):
        """Get filtered action data from engine (i.e. no zero values of lost particles)
        
        Returns
        -------
        ndarray
            filtered action array data
        """
        action = self.get_action()
        return action[action > 0]

    def get_times(self):
        """Get times from engine
        
        Returns
        -------
        ndarray
            times array
        """
        return self.times

    def get_filtered_times(self):
        """Get only loss times from engine (i.e. only loss particles)
        
        Returns
        -------
        ndarray
            filtered times array
        """
        times = self.get_times()
        return times[times < self.iterations]

    def get_survival_quota(self):
        """Get time evolution of number of survived particles
        
        Returns
        -------
        ndarray
            time evolution of survived particles
        """
        t = np.array(self.get_times())
        max_t = int(np.amax(t))
        quota = np.empty(max_t)
        for i in range(max_t):
            quota[i] = np.count_nonzero(t>i)
        return quota

    def get_lost_particles(self):
        """Get time evolution of lost particles
        
        Returns
        -------
        ndarray
            time evolution of number of lost particles
        """
        quota = self.get_survival_quota()
        return self.N - quota

    def current_binning(self, bin_size):
        """Execute current binning and computation
        
        Parameters
        ----------
        bin_size : int
            size of the binning to consider for current computation
        
        Returns
        -------
        tuple(ndarray, ndarray)
            array with corresponding sampling time (middle point), current value computed.
        """
        survival_quota = self.get_survival_quota()
        
        points = [i for i in range(0, len(survival_quota), bin_size)]
        if len(survival_quota) % bin_size == 0:
            points.append(len(survival_quota) - 1)
        t_middle = [(points[i + 1] + points[i]) *
                    0.5 for i in range(len(points) - 1)]
        currents = [(survival_quota[points[i]] - survival_quota[points[i+1]]
                     ) / bin_size for i in range(len(points) - 1)]
        return np.array(t_middle), np.array(currents)

    @staticmethod
    def generate_instance(omega_0, omega_1, omega_2, R_1, R_2, TH_MAX, barrier_radius, x_0, p_0):
       
        return symplectic_map_noise(omega_0, omega_1, omega_2, R_1, R_2, TH_MAX,  barrier_radius, x_0, p_0)

class symplectic_map_noise(symplectic_map):
    def __init__(self, omega_0, omega_1, omega_2, R_1, R_2, TH_MAX, barrier_radius, x_0, p_0):
        """Init symplectic map object!
        
        Parameters
        ----------
        object : self
            self
        omega_0 : float
            Omega 0 frequency
        omega_1 : float
            Omega 1 frequency
        omega_2 : float
            Omega 2 frequency
        R_1 : float
            inner lens radius
        R_2 : float
            outer lens radius
        TH_MAX : float
            theta max value
        barrier_radius : float
            barrier position (x coordinates!)
        x_0 : ndarray
            1D array of x initial positions
        p_0 : ndarray
            1D array of p initial values
        """
        self.omega_0 = omega_0
        self.omega_1 = omega_1
        self.omega_2 = omega_2
        self.R_1 = R_1
        self.R_2 = R_2
        self.TH_MAX = TH_MAX
        
        self.barrier_radius = barrier_radius
        self.x_0 = x_0.copy()
        self.p_0 = p_0.copy()
        self.N = len(x_0)

        self.iterations = 0

        self.x = x_0.copy()
        self.p = p_0.copy()
        self.times = np.zeros(len(self.x))

    def reset(self):
        """Reset the engine to initial conditions
        """
        self.iterations = 0
        self.x = self.x_0.copy()
        self.p = self.p_0.copy()
        self.times = np.zeros(len(self.x))

    def common_noise(self, noise_array):
        """Execute iterations with a single noise array, common for all particles.
        
        Parameters
        ----------
        noise : ndarray
            array of noise values
        """
        self.iterations += len(noise_array)
        self.x, self.p, self.times = mio.symplectic_map_common(
            self.x, self.p, self.times, noise_array, self.omega_0, self.omega_1, self.omega_2, self.R_1, self.R_2, self.TH_MAX, self.barrier_radius
        )
