from numpy import array, hstack, zeros, linspace, pi, ones
from numpy.linalg import solve
from scipy.integrate import odeint
import numpy as np
%matplotlib inline
import matplotlib.pyplot as plt
from matplotlib import cm

class Malaria_System:
    def __init__(self,pr = dict):
        #Valor dos Parametros
        self.beta_vh = pr['beta_vh'] # 0 - 1 per mosquito #
        self.kappa = pr['kappa'] # 1/11 dimensionless
        self.beta_hv = pr['beta_hv'] # Per Mosquito
        self.mu_M = pr['mu_M'] # 0.16 - 0.23 per day
        self.nu = pr['nu'] #1/15.6 +- 2.86 per day
        self.alpha = pr['alpha'] # 83 +- 48 larvae/per female mosquito
        self.mu_L = pr['mu_L'] #0.62 - 0.99 per day - Larvae mortality on ponds without fish
        self.mu_p = pr['mu_p']

        #Parametros da Vegetação

        self.r = pr['r'] # 0.5 Por Mês 
        self.gamma = pr['gamma'] # Proporção de vegetação retirada na limpeza
        self.tau = pr['tau'] # 30-60 dias
        self.H = pr['H'] # Pop. realizando Limpeza - <5%

        #Valor das Carrying Capacities Maximas
        self.K_0 = pr['K_0']
        self.K_w_max = pr['K_w_max']
        self.K_p_max = pr['K_p_max']
        self.R_0 = self.beta_hv*self.beta_vh/(self.mu_M*self.kappa)
       
    
    #Carrying capacities em função da vegetação:

    def K_p(self, v):
        #K_p is sigmoid
        #Defining Parameters 
        L = self.K_p_max # Maximum carrying capacity
        k = 5 #Steepness
        x_0 = 0.5 #Sigmoid midpoint
    
        return L/(1 + np.exp(-k*(v - x_0)))
        #return self.K_w(v)
    
    def dotK_p(self,v):
        L = self.K_p_max
        k = 5
        x_0 = 0.5
        return k*L*np.exp(-k*(v - x_0))/(1 + np.exp(-k*(v - x_0)))
    
    def K_w(self, v):
        #K_w is affine
        #Defining Parameters
        K_0 = self.K_0 # Minimum carrying capacity
        L = self.K_w_max # maximum carrying capacity
    
        return K_0 + (L - K_0)*v
    
    
    
    def _right_hand_side_v(self, x , t):
        """Returns the derivatives of the states.

        Parameters
        ----------
        x : ndarray, shape(7,1)
        The current state vector.
        t : float
        The current time.
        args : ndarray
        The constants.

        Returns
        -------
        dx : ndarray, shape(7,1)
        The derivative of the state.
    
        """

        S = x[0]
        I = x[1]
        M_S = x[2]
        M_I = x[3]
        L_p = x[4]
        L_w = x[5]
        V = x[6]
    
        S_dot = - self.beta_vh * S * M_I + self.kappa * I
        I_dot = self.beta_vh * S * M_I - self.kappa * I
        M_I_dot = self.beta_hv * I * M_S - self.mu_M * M_I
        M_S_dot = self.kappa * (L_p + L_w) - self.beta_hv * I * M_S - self.mu_M * M_S
        L_p_dot = self.alpha * (self.K_p(V)/(self.K_p(V) + self.K_w(V))) * (M_I + M_S) * (1 - L_p/self.K_p(V)) - (self.nu + self.mu_L + self.mu_p*(1 - V))*L_p
        L_w_dot = self.alpha * (self.K_w(V)/(self.K_p(V) + self.K_w(V))) * (M_I + M_S) * (1 - L_w/self.K_w(V)) - (self.nu + self.mu_L)*L_w
        V_dot = self.r*(1 - V)
    
        dx = np.array([S_dot,I_dot,M_S_dot,M_I_dot,L_p_dot,L_w_dot,V_dot])
    
        return dx

    def _jump_right_hand_side(self, x):
        """
        Returns the behaviour of the system on points of discontinuities.
    
        Parameters
        ----------
        x: ndarray, shape (7,1)
        current state vector
    
        Returns
        ----------
        dx: ndarray, shape(7,1)
        Delta between values before and after point of discontinuity
        """
        S = x[0]
        I = x[1]
        M_S = x[2]
        M_I = x[3]
        L_p = x[4]
        L_w = x[5]     
        V = x[6]
    
        S_dot = 0
        I_dot = 0
        M_I_dot = 0
        M_S_dot = 0
        L_p_dot = 0
        L_w_dot = 0
        V_dot = -self.gamma*V
    
        dx = np.array([S_dot,I_dot,M_S_dot,M_I_dot,L_p_dot,L_w_dot,V_dot])
    
        return dx
    
    def jump_odeint(self, x_initial,days_range, step):
        """
        Returns the integrated system.
        
        Parameters
        ----------
        RHS :  Right-hand side function.
    
        x_initial: vector, initial state of the system.
    
        time_vector: linspace, time over which will be integrated.
    
        jump_RHS: Discrete Delta at discontinuities.

        Returns
        -------
        y(t) : for each t in time_vector, the state y(t) of the system.
    
        """
        self.step = step
        self.time_vector =  linspace(0, days_range, int(days_range/step))
        periods = []
        for i in range(0,int(np.ceil(days_range/self.tau))):
            periods.append((i+1)*int(self.tau/step)) 

        continuous_intervals = np.split(self.time_vector,periods) #split time vector between times in jump_vector
    
        y = np.array([x_initial])
        for interval in continuous_intervals:
            if len(interval) > 0:
                y = np.append(y,odeint(self._right_hand_side_v ,y[-1],interval),axis = 0)
                y = np.append(y,np.array([y[-1] + self._jump_right_hand_side(y[-1])]),axis = 0)
        
        self.y = y[:self.time_vector.shape[0]]
        

    def _jump_right_hand_side_varying_gamma(self, x, gamma_func, t):
        """
        Returns the behaviour of the system on points of discontinuities.
    
        Parameters
        ----------
        x: ndarray, shape (7,1)
        current state vector
    
        Returns
        ----------
        dx: ndarray, shape(7,1)
        Delta between values before and after point of discontinuity
        """
        S = x[0]
        I = x[1]
        M_S = x[2]
        M_I = x[3]
        L_p = x[4]
        L_w = x[5]     
        V = x[6]
    
        S_dot = 0
        I_dot = 0
        M_I_dot = 0
        M_S_dot = 0
        L_p_dot = 0
        L_w_dot = 0
        V_dot = -gamma_func(t)*V
    
        dx = np.array([S_dot,I_dot,M_S_dot,M_I_dot,L_p_dot,L_w_dot,V_dot])
    
        return dx
    
    def jump_odeint_varying_gamma(self, x_initial,days_range, step, gamma_func):
        """
        Returns the integrated system.
        
        Parameters
        ----------
        RHS :  Right-hand side function.
    
        x_initial: vector, initial state of the system.
    
        time_vector: linspace, time over which will be integrated.
    
        jump_RHS: Discrete Delta at discontinuities.

        Returns
        -------
        y(t) : for each t in time_vector, the state y(t) of the system.
    
        """
        self.step = step
        self.time_vector =  linspace(0, days_range, int(days_range/step))
        periods = []
        for i in range(0,int(np.ceil(days_range/self.tau))):
            periods.append((i+1)*int(self.tau/step)) 

        continuous_intervals = np.split(self.time_vector,periods) #split time vector between times in jump_vector
    
        y = np.array([x_initial])
        for interval in continuous_intervals:
            if len(interval) > 0:
                y = np.append(y,odeint(self._right_hand_side_v ,y[-1],interval),axis = 0)
                y = np.append(y,np.array([y[-1] + self._jump_right_hand_side_varying_gamma(y[-1],gamma_func,interval[-1])]),axis = 0)
        
        self.y = y[:self.time_vector.shape[0]]


        
    def Veg(self, t):
        return 1 - (self.gamma * np.exp(-self.r* (t % self.tau)))/(1 - (1 - self.gamma)*np.exp(-self.r*self.tau))
    
    def Basic_offspring_num(self):
        return (self.alpha*self.nu)/((self.nu + self.mu_L)*self.mu_M)
    
    def Upp_thresh_Kw(self):
        t = linspace(0, self.tau, self.tau*10000)
        K_0 = self.K_0 # Minimum carrying capacity
        L = self.K_w_max # maximum carrying capacity
        dotK = (L - K_0)*self.r*(1 - self.Veg(t))/(self.K_w(self.Veg(t)))
        
        return (np.max(dotK) + self.nu + self.mu_L)*(np.max(self.K_w(self.Veg(t)) + self.K_p(self.Veg(t))))/np.min(self.K_w(self.Veg(t)))

    def Upp_thresh_Kp(self):
        t = linspace(0, self.tau, self.tau*10000)
        
        K_0 = self.K_0 # Minimum carrying capacity
        L = self.K_w_max # maximum carrying capacity

        dot_term = self.dotK_p(self.Veg(t))*self.r*(1 - self.Veg(t))/(self.K_p(self.Veg(t)))
        
        return (np.max(dot_term) + self.nu + self.mu_L + self.mu_p*(1 - np.min(self.Veg(t))))*(np.max(self.K_w(self.Veg(t)) + self.K_p(self.Veg(t))))/np.min(self.K_p(self.Veg(t)))

    def Low_thresh_Kw(self):
        t = linspace(0, self.tau, self.tau*10000)
        
        K_0 = self.K_0 # Minimum carrying capacity
        L = self.K_w_max # maximum carrying capacity
        dotK = (L - K_0)*self.r*(1 - self.Veg(t))/(self.K_w(self.Veg(t)))
        
        return (np.min(dotK) + self.nu + self.mu_L)*(np.min(self.K_w(self.Veg(t)) + self.K_p(self.Veg(t))))/np.max(self.K_w(self.Veg(t)))

    def Low_thresh_Kp(self):
        
        t = linspace(0, self.tau, self.tau*10000)
        
        K_0 = self.K_0 # Minimum carrying capacity
        L = self.K_w_max # maximum carrying capacity

        dot_term = self.dotK_p(self.Veg(t))*self.r*(1 - self.Veg(t))/(self.K_p(self.Veg(t)))
        
        return (np.min(dot_term) + self.nu + self.mu_L + self.mu_p*(1 - np.max(self.Veg(t))))*(np.min(self.K_w(self.Veg(t)) + self.K_p(self.Veg(t))))/np.max(self.K_p(self.Veg(t)))


    def Mosquito_threshold(self):        
        
        
        thres = self.Basic_offspring_num()
        
        Upp_Kp = self.Upp_thresh_Kp()
        Upp_Kw = self.Upp_thresh_Kw()
        
        Low_Kp = self.Low_thresh_Kp()
        Low_Kw = self.Low_thresh_Kw()

        print("Basic Offspring Number = ",thres, "\n \n Upper Bound for Kp,Kw = ",(Upp_Kp,Upp_Kw), "\n\n Lower Bound Kp,Kw = ", (Low_Kp,Low_Kw))
        
        A = Upp_Kp >= thres and Upp_Kw >= thres
        B = Low_Kp < thres or Low_Kw < thres 
        if A: print("Limit to zero");   return False
        if B: print("Limit to periodic solution");  return True
        if not A and not B: print("Inconclusive Result"); return None


    def dis_min_thresh(self):
        bigM = self.y.T[2,-int(self.tau/self.step):] + self.y.T[3,-int(self.tau/self.step):] 
        dotM = self.kappa*(self.y.T[4,-int(self.tau/self.step):] + self.y.T[5,-int(self.tau/self.step):])/bigM - self.mu_M
        return self.beta_hv*self.beta_vh*np.min(bigM)/(self.kappa*(self.mu_M - np.max(dotM)))
        
    def dis_max_thresh(self):
        bigM = self.y.T[2,-int(self.tau/self.step):] + self.y.T[3,-int(self.tau/self.step):] 
        dotM = self.kappa*(self.y.T[4,-int(self.tau/self.step):] + self.y.T[5,-int(self.tau/self.step):])/bigM - self.mu_M
        return self.beta_hv*self.beta_vh*np.max(bigM)/(self.kappa*(self.mu_M - np.min(dotM)))
        
    def Disease_threshold(self):
        min_thresh = self.dis_min_thresh()
        max_thresh = self.dis_max_thresh()
        print("min = ", min_thresh,"\n")
        print("max = ", max_thresh,"\n")
        if min_thresh > 1:
            print("periodic solution")
        if max_thresh <= 1:
            print("Solution converges to zero")
        if min_thresh <= 1 and max_thresh > 1:
            print("Inconclusive result")
        return
      
    def R_0_t(self):
        bigM = self.y.T[2] + self.y.T[3]
        return self.beta_hv*self.beta_vh*bigM/(self.kappa*self.mu_M)
    
    def Plot(self):
        
        t = self.time_vector
        fig, axs = plt.subplots(4,1, figsize=(21, 30))
        
        font_size = 30
        tick_label_size = 25
        
        axs[0].plot(t, self.y[:,6], label='$V(t)$')
        axs[0].plot(t, self.y[:,4], label='$L_p(t)$')
        axs[0].plot(t, self.y[:,5], label = '$L_w(t)$')
        axs[0].set_xlabel('time (days)', fontsize = font_size)
        axs[0].set_ylabel('', fontsize = font_size)
        axs[0].grid(True)
        axs[0].tick_params(axis = 'both', labelsize = tick_label_size)
        axs[0].margins(x = 0)
        axs[0].legend(fontsize = font_size, loc='upper right')


        axs[1].plot(t, self.y[:,2], label = '$M_S(t)$')
        axs[1].plot(t, self.y[:,3], label = '$M_I(t)$')
        axs[1].set_xlabel('time (days)', fontsize = font_size)
        axs[1].set_ylabel('', fontsize = font_size)
        axs[1].grid(True)
        axs[1].tick_params(axis = 'both', labelsize = tick_label_size)
        axs[1].margins(x = 0)
        axs[1].legend(fontsize = font_size, loc='upper right')

        axs[2].plot(t, self.y[:,0], label = '$S(t)$')
        axs[2].plot(t, self.y[:,1], label = '$I(t)$')
        axs[2].set_xlabel('time (days)', fontsize = font_size)
        axs[2].set_ylabel('', fontsize = font_size)
        axs[2].grid(True)
        axs[2].tick_params(axis = 'both', labelsize = tick_label_size)
        axs[2].margins(x = 0)
        axs[2].legend(fontsize = font_size, loc='upper right')

        axs[3].plot(t, self.R_0_t(), label = '$R(t)$')
        axs[3].set_xlabel('time (days)', fontsize = font_size)
        axs[3].set_ylabel('', fontsize = font_size)
        axs[3].grid(True)
        axs[3].tick_params(axis = 'both', labelsize = tick_label_size)
        axs[3].margins(x = 0)
        axs[3].legend(fontsize = font_size, loc='upper right')

        plt.show()
        return

    def Plot_column(self, axs, column):
        
        t = self.time_vector

        font_size = 30
        tick_label_size = 25

        axs[0, column].plot(t, self.y[:,4], label='$L_p(t)$', dashes = [6,2])
        axs[0, column].plot(t, self.y[:,5], label = '$L_w(t)$')
        axs[0, column].set_xlabel('', fontsize = font_size)
        axs[0, column].set_ylabel('', fontsize = font_size)
        axs[0, column].grid(True)
        axs[0, column].tick_params(axis = 'both', labelsize = tick_label_size)
        axs[0, column].margins(x = 0)
        #axs[0, column].legend(fontsize = font_size, loc='upper right')


        axs[1, column].plot(t, self.y[:,2], label = '$M_S(t)$')
        axs[1, column].plot(t, self.y[:,3], label = '$M_I(t)$', dashes = [6,2])
        axs[1, column].set_xlabel('', fontsize = font_size)
        axs[1, column].set_ylabel('', fontsize = font_size)
        axs[1, column].grid(True)
        axs[1, column].tick_params(axis = 'both', labelsize = tick_label_size)
        axs[1, column].margins(x = 0)
        #axs[1, column].legend(fontsize = font_size, loc='upper right')

        axs[2, column].plot(t, self.y[:,0], label = '$S(t)$')
        axs[2, column].plot(t, self.y[:,1], label = '$I(t)$', dashes = [6,2])
        axs[2, column].set_xlabel('', fontsize = font_size)
        axs[2, column].set_ylabel('', fontsize = font_size)
        axs[2, column].grid(True)
        axs[2, column].tick_params(axis = 'both', labelsize = tick_label_size)
        axs[2, column].margins(x = 0)
        #axs[2, column].legend(fontsize = font_size, loc='upper right')

        #axs[3, column].plot(t,np.ones(len(t)),color = 'r')
        axs[3, column].plot(t, self.y[:,6], label='$V(t)$',dashes = [6,2])
        axs[3, column].plot(t, self.R_0_t(), label = '$R(t)$')
        axs[3, column].set_xlabel('', fontsize = font_size)
        axs[3, column].set_ylabel('', fontsize = font_size)
        axs[3, column].grid(True)
        axs[3, column].set_yticks([1,5,10,15])
        axs[3, column].tick_params(axis = 'both', labelsize = tick_label_size)
        axs[3, column].margins(x = 0)
        #axs[3, column].legend(fontsize = font_size, loc='upper right')

        return

    def plot_row(self, axs, row):
        
        t = self.time_vector

        font_size = 30
        tick_label_size = 25

        axs[row, 0].plot(t, self.y[:,4], label='$L_p(t)$', dashes = [6,2])
        axs[row, 0].plot(t, self.y[:,5], label = '$L_w(t)$')
        axs[row, 0].set_xlabel('', fontsize = font_size)
        axs[row, 0].set_ylabel('', fontsize = font_size)
        axs[row, 0].grid(True)
        axs[row, 0].tick_params(axis = 'both', labelsize = tick_label_size)
        axs[row, 0].margins(x = 0)
        #axs[row, 0].legend(fontsize = font_size, loc='upper right')


        axs[row, 1].plot(t, self.y[:,2], label = '$M_S(t)$')
        axs[row, 1].plot(t, self.y[:,3], label = '$M_I(t)$', dashes = [6,2])
        axs[row, 1].set_xlabel('', fontsize = font_size)
        axs[row, 1].set_ylabel('', fontsize = font_size)
        axs[row, 1].grid(True)
        axs[row, 1].tick_params(axis = 'both', labelsize = tick_label_size)
        axs[row, 1].margins(x = 0)
        #axs[row, 1].legend(fontsize = font_size, loc='upper right')

        axs[row, 2].plot(t, self.y[:,0], label = '$S(t)$')
        axs[row, 2].plot(t, self.y[:,1], label = '$I(t)$', dashes = [6,2])
        axs[row, 2].set_xlabel('', fontsize = font_size)
        axs[row, 2].set_ylabel('', fontsize = font_size)
        axs[row, 2].grid(True)
        axs[row, 2].tick_params(axis = 'both', labelsize = tick_label_size)
        axs[row, 2].margins(x = 0)
        #axs[row, 2].legend(fontsize = font_size, loc='upper right')

        #axs[row, 3].plot(t,np.ones(len(t)),color = 'r')
        axs[row, 3].plot(t, self.y[:,6], label='$V(t)$',dashes = [6,2])
        axs[row, 3].plot(t, self.R_0_t(), label = '$R(t)$')
        axs[row, 3].set_xlabel('', fontsize = font_size)
        axs[row, 3].set_ylabel('', fontsize = font_size)
        axs[row, 3].grid(True)
        axs[row, 3].set_yticks([1,5,10,15])
        axs[row, 3].tick_params(axis = 'both', labelsize = tick_label_size)
        axs[row, 3].margins(x = 0)
        #axs[row, 3].legend(fontsize = font_size, loc='upper right')

        return