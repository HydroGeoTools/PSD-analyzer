import numpy as np
import scipy.optimize

# d'après:
# https://publications.polymtl.ca/3943/1/2019_RoselineTaillonLevesque.pdf


g = 9.81 # gravity (m/s^2)
rho_w = 998 #water density (20°C)
gamma_w = rho_w * g # unit weight of water (N/m^3)
vu_w = 1.005e-3 #kinematic viscosity (N s /m^2 @ 20°C)
mu_w = vu_w * rho_w #dynamic viscosity
nu_w = g * rho_w # volumic weight 
sigma_w = 0.07275 # water surface tension (N/m @ 20°C)


###
# Curve fitting de la PSD
###

def Fredlund_PSD(d, log_a_gr, log_n_gr, log_m_gr, log_dr, log_dm):
    #https://cdnsciencepub.com/doi/10.1139/t00-015
    a_gr = 10**log_a_gr
    n_gr = 10**log_n_gr
    m_gr = 10**log_m_gr
    dr = 10**(-log_dr)
    dm = 10**(-log_dm)
    res = 1 / np.log(np.exp(1)+(a_gr/d)**n_gr)**m_gr
    res *= 1 - (np.log(1+dr/d)/np.log(1+dr/dm))**7
    return 100*res

def MSE(x, model, xdata, ydata):
    y_th = model(xdata, *x)
    residual = y_th - ydata
    MSE = 1/len(xdata) * np.sum(residual**2)
    return MSE

def R2(MSE, ydata):
    y_mean = np.mean(ydata)
    SStot = np.sum((ydata - y_mean)**2)
    SSres = MSE * len(ydata)
    R2 = 1 - SSres / SStot
    return R2


def Dx(val, func, params):
    """
    Return the D_x particle size such as D_10, D_30, ...
    """
    x0 = -2.
    res = scipy.optimize.root(lambda n: func(10**n, *params) - val, x0)
    print(val, res)
    return float(10**res.x)


def fit(xdata, ydata):
    #from https://cdnsciencepub.com/doi/10.1139/t00-015, figure 11
    func = Fredlund_PSD
    bounds=[
        (-6,1), #log_a_gr_bounds
        (0,20), #log_n_gr_bounds
        (-1,5), #log_m_gr_bounds
        (-5,5), #log_dr_bounds
        (0,6), #log_dm_bounds
    ]
    res = scipy.optimize.dual_annealing(MSE, bounds=bounds, args=[func, xdata, ydata], maxiter=1000, initial_temp=1e4, restart_temp_ratio=1e-4)
    return res

def surface_specific(rho_s, func):
    # return specific surface based on the PSD
    Ss = 6/rho_s * scipy.integrate.quad()
    return Ss


###
# Prediction of saturated hydraulic conductivity
###  

# See also: https://link.springer.com/article/10.1007/s10064-012-0418-7

def modified_Kozeny_Carman_granular(e, C_U, D_10):
    # https://link.springer.com/article/10.1023/A:1016046214724
    # note: D_10 should be in m!
    x = 2
    C_G = 0.1 
    K = C_G * gamma_w / mu_w * e**(3+x) / (1+e) * C_U**(1/3) * D_10**2 #equation 13
    return K

def modified_Kozeny_Carman_plastic(e, rho_s, w_l):
    # https://link.springer.com/article/10.1023/A:1016046214724
    x = np.min([2,7.7 * w_l**(-0.15) - 3]) #equation 6
    C_P = 5.6 #g^2 / m^4
    chi = 1.5
    K = C_P * gamma_w / mu_w * e**(3+x) / (1+e) * 1 / (e*rho_s**2 * w_l**(2*chi)) #equation 17
    return K


###
# Prediction of the WRC
###  

def modified_Kovacs(phi, e, D_60, D_10):
    # https://cdnsciencepub.com/doi/10.1139/t03-054
    a_c = 0.1
    C_U = D_60 / D_10
    #for low-plasticity low-cohesive soil
    alpha = 10
    beta_w = 0
    b = alpha * sigma_w * np.cos(beta) / ((1.17*np.log10(C_U)+1) * gamma_w)
    h_co = b / (e * D_10)
    
    D_H = (1 + 1.17 * np.log10(C_U)) * D_10
    phi_r = 0.42 / (e * D_H)**1.26
    m = 1/C_U
    S_c = 1-((h_co/phi)**2+1)**m * np.exp(-m*(h_co/phi)**2)
    C_phi = 1 - np.log(1 + phi/phi_r) / np.log(1+phi_0/phi_r)
    S_a = a_c * C_phi * (h_co/phi_n)**(2/3) / (e**(1/3) * (phi/phi_n)**(1/6))
    S_r = 1 - np.where(1-Sa <= 0, 0, 1-Sa) * (1-Sc) #relative saturation
    return S_r
