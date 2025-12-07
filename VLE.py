import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

import scipy as sc

"""
This module contains functions and constants to model a batch distillation of an ethanol-water mixture.

Numpy is used for numerical operations, matplotlib for plotting, pandas for data manipulation,
and scipy for scientific computations. Pls update the path to the excel files on your computer and
make sure you have the required libraries installed.

pip install numpy matplotlib pandas scipy : is a command to install the required libraries if not already installed.
Don't forget to run the code in an environment that supports these libraries.
https://numpy.org/
https://scipy.org/
https://matplotlib.org/
https://pandas.pydata.org/


"""






""" 
Functions 
"""



def L_rg(x_eth,x_eth_init,y_eth_interp,L_init):#OK 
    """_summary_
        This function gives the numerical solution of the Rayleigh equation
        Assuming real gas with y_eth_interp(x_eth).
        Must be computed in mole as we derived the equation in molar fraction.
    Args:
        x_eth (_type_): _numpy-vector/scalar_
        x_eth_init (_type_, optional):  Defaults to x_eth_init.
        y_eth_rg (_type_, optional):  Defaults to y_eth_rg.
        L_init (_type_, optional): Defaults to L0.

    Returns:
        _type_: _numpy-vector/scalar_
    """
    ln_L = np.log(L_init) + sc.integrate.quad(lambda x : 1/(y_eth_interp(x)-x),x_eth_init,x_eth)[0]
    return np.exp(ln_L)

def vapor_mass_fraction_to_distillate_mass_fraction(Distillate,w_eth_vapor_interp_mass_D,D_init = 0.0):#OK
        """_summary_
            This function gives the mass fraction of ethanol in the distillate
            as a function of the mass of distillate collected.
            It is obtained by integrating the mass fraction of ethanol in the vapor phase
            over the mass of distillate collected.

        Args:
            Distillate (_type_): _numpy-vector_

        Returns:
            _type_: __numpy-vector_
        """
        if Distillate == 0.0:
            Distillate += 0.000001
            
        return sc.integrate.quad(lambda d : w_eth_vapor_interp_mass_D(d),D_init,Distillate)[0]/Distillate

def h_vap(h_vap_w, h_vap_eth, w_liq, T):
    """This function aims to compute the specific enthalpy of vaporization of our distillate vapor assuming ideal gas.

    Args:
        h_vap_w (func): specific enthalpy of vaporization of water
        h_vap_eth (func): specific enthalpy of vaporization of ethanol
        w_liq (numpy 1d array): array of the weight percentage of ethanol in liquid phase
        T (numpy 1d array): array of the temperature

    Returns:
        numpy 1d array : specific enthalpy of vaporization 
    """
    return w_liq*h_vap_eth(T) + (1-w_liq)*h_vap_w(T)

def dT_dD (T,D,x):
    """This function aims to compute the numerical derivative of dT/dD

    Args:
        T (numpy 1d array): array of temperatures
        D (numpy 1d array): array of distillate mass
        x (numpy 1d array): array of composition of liquid 

    Returns:
        _type_: _description_
    """
    return np.gradient(T,x)/np.gradient(D,x)

def m_dot( Q_dot , h_vap , Cp, dT_dD, L_0, D, T):
    """_summary_
            This function aims to compute the mass flow at each step of our distillation based 
            on an energy balance on a signle batch distillation.
            
        Args:
            Q_dot : a constant depending on the heater [kW]
            h_vap : enthalpy of vaporization of our mixture, is a function of the temperature and 
                    the composition of the liquid and vapor phase. [kJ/kg]
            Cp    : The heat capacity of our different liquid phase (assumed constant) [kJ/kg*K]
            T     : Temperature of the liquid phase as a function of the distillate mass [K]
            L_0   : Initial mass of the liquid phase [kg]
            D     : Mass of our distillate [kg]

        Returns:
            m_dot : mass flow [kg/s]
    
    """
    
    
    return Q_dot/(h_vap - Cp*(273.15+T) + (L_0 - D)*Cp*dT_dD)


def time_of_distillation(m_dot,D,D_init = 0.0):
    """_summary_
            This function aims to compute the time needed to perform the distillation by solving the differential equation :
            dD/dt = m_dot --> t = int(1/m_dot dD)
            
        Args:
            m_dot : mass flow rate [kg/s]
            D     : mass of our distillate [kg]


        Returns:
            m_dot : mass flow [kg/s]
    """
    m_dot_interp = sc.interpolate.interp1d(D, m_dot, kind='cubic')

    return sc.integrate.quad(lambda d: 1/m_dot_interp(d),D_init,D[-1])[0]



"""
Table of conversion

"""
#PLS UPDATE THE PATH TO THE EXCEl FILE ON YOUR COMPUTER
conv_table = pd.read_excel(r"C:\Users\guylo\Documents\Unif\Bac 3\Q1\Bachelor\Distillation_team\Table desnity - Water-Sugar & Water Ethanol.xlsx",sheet_name="Water Ethanol")
volume_percentage = conv_table['%_Veth'].to_numpy()
rho_mixt = conv_table['Density '].to_numpy()
mole_percentage = conv_table['%_neth'].to_numpy()
mass_percentage = conv_table['%_meth'].to_numpy()

#Creation of interpolation fonctions to perform conversions

rho_mixt_interp_v = sc.interpolate.interp1d(volume_percentage,rho_mixt,kind = 'cubic')#kg/m³
rho_mixt_interp_n = sc.interpolate.interp1d(mole_percentage,rho_mixt,kind = 'cubic')#kg/m³
rho_mixt_interp_m = sc.interpolate.interp1d(mass_percentage,rho_mixt,kind = 'cubic')#kg/m³


mole_percentage_interp_v = sc.interpolate.interp1d(volume_percentage,mole_percentage,kind = 'cubic')#%
mole_percentage_interp_m =sc.interpolate.interp1d(mass_percentage,mole_percentage,kind = 'cubic')#%

mass_percentage_interp_v = sc.interpolate.interp1d(volume_percentage,mass_percentage,kind = 'cubic')#%
mass_percentage_interp_n = sc.interpolate.interp1d(mole_percentage,mass_percentage,kind = 'cubic')#%

volume_percentage_interp_n = sc.interpolate.interp1d(mole_percentage,volume_percentage,kind = 'cubic')
volume_percentage_interp_m = sc.interpolate.interp1d(mass_percentage,volume_percentage,kind = 'cubic')

"""

Defininition of the system constants and initial conditions : IC


"""
volumique_percentage_init = 15#%


MM_w = 18.02 #[kg/kmol]
MM_eth = 46.07 #[kg/kmol]
rho_eth = 789 #[kg/m³]

L0 = 0.200449143




x_eth_init = 0.183884096




V_eth_init =2.151416445/1000 #m³







#THERMODYNAMICS VARIABLE FROM EES

Cp_w = 4.200 #[kJ/kg] assumed constant because Cp = [4.18,4.214] [25,100]
Cp_eth = 2.917 #[kJ/kg] assumed constant in our range bc Cp = [2.434,2.917] [25,100] ATT en phase liq

df_h_vap = pd.read_excel(r"C:\Users\guylo\Documents\Unif\Bac 3\Q1\Bachelor\Distillation_team\specific_enthalpy.xlsx",sheet_name="specific enthalpy of vap")
df_h_vap_w = df_h_vap['h_vap_w'].to_numpy()#[kJ/kg]
df_h_vap_eth = df_h_vap['h_vap_eth'].to_numpy()#[kJ/kg]
df_h_vap_T = df_h_vap['Temperature'].to_numpy()
h_vap_eth = sc.interpolate.interp1d(df_h_vap_T,df_h_vap_eth,kind = 'cubic')
h_vap_w = sc.interpolate.interp1d(df_h_vap_T,df_h_vap_w,kind = 'cubic')

df_rg = pd.read_excel(r"C:\Users\guylo\Documents\Unif\Bac 3\Q1\Bachelor\Distillation_team\VLE_x_y.xlsx",sheet_name="data from literature")
df_x_1 = df_rg['x_1'].to_numpy()
df_y_1 = df_rg['y_1'].to_numpy()
df_T = df_rg['T_B in °C'].to_numpy() 
df_aplha = df_rg['alpha_1,2'].to_numpy()

Temperature_x = sc.interpolate.interp1d(df_x_1, df_T, kind='cubic') # function of x
Temperature_y = sc.interpolate.interp1d(df_y_1, df_T, kind='cubic') # function of y

y_eth = sc.interpolate.interp1d(df_x_1, df_y_1, kind='cubic') # function of x



"""
COMPUTING : MASS EQUILIBRIUM (where to stop distillation)
"""



"Here is the different percentages of ethanol in liquid phase"
x1 = np.linspace(x_eth_init, 0.039, 500)
x2 = np.logspace(np.log10(0.04), -16, 500)

x_eth = np.concatenate([x1, x2]) 

vol_x_eth = volume_percentage_interp_n(x_eth*100)#%
mass_x_eth = mass_percentage_interp_n(x_eth*100)#%
print(mass_x_eth[0])

"Here is the different percentages of ethanol in the vapor phase and so in the distillate"
y_eth_init = y_eth(x_eth_init)
y_eth_val = y_eth(x_eth)
vol_y_eth = volume_percentage_interp_n(y_eth_val*100)#%
mass_y_eth = mass_percentage_interp_n(y_eth_val*100)#%

"Rayleigh equation solution in molar; volume and mass"
L_rg_val = np.zeros_like(x_eth)
for i,x in enumerate(x_eth) :
    L_rg_val[i] = L_rg(x,x_eth_init,y_eth,L0) #kmol

vol_L_rg_val = L_rg_val*(MM_eth*mole_percentage_interp_v(vol_x_eth)/100 + MM_w*(1-mole_percentage_interp_v(vol_x_eth)/100))/rho_mixt_interp_v(vol_x_eth) 
mass_L_rg_val = L_rg_val*(MM_eth*mole_percentage_interp_v(vol_x_eth)/100 + MM_w*(1-mole_percentage_interp_v(vol_x_eth)/100)) #kg

"Amount of distillate colected in molar; volume and mass"
D = L0 - L_rg_val #kmol
vol_D = vol_L_rg_val[0] - vol_L_rg_val
mass_D = mass_L_rg_val[0] - mass_L_rg_val
print(mass_D)









"Computing the cumulative amount of ethanol in the distillate collected"
mass_fraction_in_vapor_interp_D = sc.interpolate.interp1d(mass_D,mass_y_eth,kind = 'cubic')
mass_fraction_distillate = np.zeros_like(mass_y_eth) # Will be used to determine the cumulative mass fraction in distillate

for i,d in enumerate(mass_D):
    mass_fraction_distillate[i] = vapor_mass_fraction_to_distillate_mass_fraction(d,mass_fraction_in_vapor_interp_D)
    
vol_percentage_distillate = volume_percentage_interp_m(mass_fraction_distillate)
mole_percentage_distillate = mole_percentage_interp_m(mass_fraction_distillate)
vol_eth = vol_percentage_distillate*vol_D/100
rho_D = rho_mixt_interp_v(vol_percentage_distillate)
rho_L = rho_mixt_interp_v(vol_x_eth)
"""

COMPUTING : ENERGY BALANCE (when to stop distillation) 
We assume ideal gas for the vapor phase enthalpy calculation
We assume constant Cp for the liquid phase enthalpy calculation
We assume constant Q_dot provided by the heater
We assume no heat loss to the environment

"""
Q_dot= 1.4 #kW
T = Temperature_x(x_eth)  #C
h_vap_mixt = h_vap(h_vap_w, h_vap_eth, mass_percentage_interp_v(x_eth*100)/100, T) #kJ/kg
dT_dD_vals = dT_dD (T,mass_D,x_eth) #C/kg as glucose don't evaporate, the derivative is the same as for water-ethanol
mass_tripartite_init = V_batch*rho_real #kg
Cp = Cp_eth*(V_eth_init*rho_eth/mass_tripartite_init) + Cp_w*(1-V_eth_init*rho_eth/mass_tripartite_init) #kJ/kg*C
""" 
It seems that Cp of water + sugar is less than Cp of pure water so we will assume Cp water to be conservative

"""
mass_tripartite_init = V_batch*rho_real #kg
m_dot_vals = m_dot( Q_dot , h_vap_mixt , Cp, dT_dD_vals, mass_tripartite_init, mass_D, T) #kg/s
vol_m_dot_vals = m_dot_vals/(rho_mixt_interp_v(vol_x_eth)) #m³/s
time_of_distillation_needed = time_of_distillation(m_dot_vals,mass_D,D_init=mass_D[0]) #s
print("Time needed for the distillation : ",time_of_distillation_needed/60," minutes")


if __name__ == "__main__" :
  
    name = "all"
    
    "volatility_curve\n"
    "dew_point_bubble_point\n"
    "x_y_curve\n"
    "x_vs_volume\n"
    "y_vs_volume\n"
    "vol_x_vs_volume\n"
    "vol_y_vs_volume\n"
    "cumulative_vol_y_vs_volume\n"
    "m_dot_vs_volume\n"
    "time_needed\n"
    "v_eth/v_eth_init\n"
    "all"
    
    if name == "volatility_curve" or name == "all":
       
        figure, ax = plt.subplots()
        x_eth_prime = df_x_1[1:-1]
        alpha = df_aplha[1:-1]
        ax.plot(x_eth_prime, alpha, label='Volatility curve', color='blue')

        ax.set_xlabel(r'$x_{eth}$ (%mol)')
        ax.set_ylabel(r'$\alpha$')
        ax.set_title('Volatility curve of ethanol in water at 1 atm')
        ax.grid()
        ax.legend()
        plt.show()
    
    if name == "dew_point_bubble_point" or name == "all":
        figure, ax = plt.subplots()
        ax.plot(df_x_1, Temperature_x(df_x_1), label='Bubble point', color='blue')
        ax.plot(df_y_1, Temperature_y(df_y_1), label='Dew point', color='orange')
        ax.scatter(x_eth_init, Temperature_x(x_eth_init), label='x0  ', color='red')
        ax.scatter(y_eth_init, Temperature_y(y_eth_init), label='y0', color='red')
        ax.set_xlabel('molar fraction of ethanol')
        ax.set_ylabel('Temperature (°C)')
        ax.set_title('Bubble and dew point curves of ethanol in water at 1 atm')
        ax.grid()
        ax.legend()
        plt.show()
        
    if name == "x_y_curve" or name == "all":
        figure, ax = plt.subplots()
        ax.plot(df_x_1, df_y_1, label='x-y curve', color='blue')
        ax.plot([0,1],[0,1], label='y=x', color='orange', linestyle='--')
        ax.scatter(x_eth_init, y_eth_init, label='initial point', color='red')
        ax.set_xlabel('x_eth (%mol)')
        ax.set_ylabel('y_eth (%mol)')
        ax.set_title('x-y curve of ethanol in water at 1 atm')
        ax.grid()
        ax.legend()
        plt.show()
        
    
        figure,ax = plt.subplots()
        ax.plot(mass_D,1/m_dot_vals, label=r'$\frac{1}{\dot{m}}$', color='red') 
        ax.set_xlabel('mass of distillate (kg)')
        ax.set_ylabel(r'$\frac{1}{\dot{m}}$ (s/kg)')
        ax.set_title(r'Plot of $\frac{1}{\dot{m}}$ as a function of mass of distillate collected')
        x_fill = mass_D
        y_fill = 1/m_dot_vals
        plt.fill_between(x_fill, y_fill, color='skyblue', alpha=0.4, label=fr'time of distillation : {time_of_distillation_needed /60:.2f} min')

        ax.grid()
        ax.legend()
        plt.show()
    
    if name == "x_vs_volume" or name == "all":
        figure,ax = plt.subplots()
        ax.plot(vol_L_rg_val*1000,x_eth*100, label = "Liquid phase", color = 'blue')
        ax.plot(vol_D*1000,x_eth*100, label = 'Distillate', color = 'orange')
        ax.set_xlabel('volume (L)')
        ax.set_ylabel(r' $x_{eth}(\%)$')
        
        ax.grid()
        ax.legend()
        plt.show()
    
    if name == "y_vs_volume" or name == "all":
        figure,ax = plt.subplots()
        ax.plot(vol_L_rg_val*1000,y_eth_val*100, label = "Liquid phase", color = 'blue')
        ax.plot(vol_D*1000,y_eth_val*100, label = 'Distillate', color = 'orange')
        ax.set_xlabel('volume (L)')
        ax.set_ylabel(r' $y_{eth}(\%)$')
        
        ax.grid()
        ax.legend()
        plt.show()
        
    if name == 'vol_x_vs_volume' or name == 'all':
        figure,ax = plt.subplots()
        ax.plot(vol_L_rg_val*1000,vol_x_eth, label = "Liquid phase", color = 'blue')
        ax.plot(vol_D*1000,vol_x_eth, label = 'Distillate', color = 'orange')
        ax.set_xlabel('volume (L)')
        ax.set_ylabel(r' Volume percentage in the liquid phase(%)')
        
        ax.grid()
        ax.legend()
        plt.show()
        
    if name == 'vol_y_vs_volume' or name == 'all':
        figure,ax = plt.subplots()
        ax.plot(vol_L_rg_val*1000,vol_y_eth, label = "Liquid phase", color = 'blue')
        ax.plot(vol_D*1000,vol_y_eth, label = 'Distillate', color = 'orange')
        ax.set_xlabel('volume (L)')
        ax.set_ylabel(r' Volume percentage in the distillate(%)')
        
        ax.grid()
        ax.legend()
        plt.show()
        
    if name == 'cumulative_vol_y_vs_volume' or name == 'all':
        figure,ax = plt.subplots()
        ax.plot(vol_D*1000,vol_percentage_distillate, label = 'Distillate', color = 'orange')
        ax.plot(vol_L_rg_val*1000,vol_percentage_distillate,label = "Liquid phase", color = 'blue')
        ax.plot(vol_D*1000,volumique_percentage_init*np.ones_like(vol_D), label = 'Initial volume percentage : 15%', color = 'red', linestyle='--')
        ax.set_xlabel('volume of distillate collected (L)')
        ax.set_ylabel(r' Cumulative volume percentage in the distillate(%)')
        
        ax.grid()
        ax.legend()
        plt.show()
        
    if name == 'm_dot_vs_volume' or name == 'all':
        figure,ax = plt.subplots()
        ax.plot(vol_D*1000,m_dot_vals, label = 'Distillate', color = 'orange')
        ax.set_xlabel('volume of distillate collected (L)')
        ax.set_ylabel(r' $\dot{m} (kg/s)$')
        
        ax.grid()
        ax.legend()
        plt.show()
        
    if name == 'time_needed' or name == 'all':
        print("Time needed for the distillation : ",time_of_distillation_needed/60," minutes")
        figure, ax = plt.subplots()
        ax.plot(mass_D,1/m_dot_vals, label=r'$\frac{1}{\dot{m}}$', color='red') 
        ax.set_xlabel('mass of distillate (kg)')
        ax.set_ylabel(r'$\frac{1}{\dot{m}}$ (s/kg)')
        ax.set_title(r'Plot of $\frac{1}{\dot{m}}$ as a function of mass of distillate collected')
        x_fill = mass_D
        y_fill = 1/m_dot_vals
        plt.fill_between(x_fill, y_fill, color='skyblue', alpha=0.4, label=fr'time of distillation : {time_of_distillation_needed /60:.2f} min')

        ax.grid()
        ax.legend()
        plt.show()
        
    if name == 'v_eth/v_eth_init' or name == 'all':
        figure,ax = plt.subplots()
        ax.plot(vol_D*1000,vol_eth/ (V_eth_init), label = 'Distillate', color = 'orange')
        ax.set_xlabel('volume of distillate collected (L)')
        ax.set_ylabel(r' $\frac{V_{ethanol}}{V_{ethanol,init}}$')
        idx = np.where( vol_percentage_distillate <= 45.0)[0]
        ax.axhline( y = vol_eth[idx][0]/ (V_eth_init), color = 'red', linestyle='--', label='45% vol')
        ax.axvline( x = vol_D[idx][0]*1000, color = 'red', linestyle='--', label=f'Volume needed for 45% vol : {vol_D[idx][0]*1000:.2f} L')
        ax.grid()
        ax.legend()
        plt.show()
    else :
        print(
                "Pls write correctly :\n"
                "volatility_curve\n"
                "dew_point_bubble_point\n"
                "x_y_curve\n"
                "x_vs_volume\n"
                "y_vs_volume\n"
                "vol_x_vs_volume\n"
                "vol_y_vs_volume\n"
                "cumulative_vol_y_vs_volume\n"
                "m_dot_vs_volume\n"
                "time_needed\n"
                "v_eth/v_eth_init\n"
                "all")

    
