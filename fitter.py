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

model_output_name = ["a_gr", "n_gr", "m_gr", "dr", "dm"]


###
# Curve fitting de la PSD
###

def Fredlund_PSD(d, log_a_gr, log_n_gr, log_m_gr, log_dr, log_dm):
    #https://cdnsciencepub.com/doi/10.1139/t00-015
    a_gr = 10**log_a_gr
    n_gr = 10**log_n_gr
    m_gr = 10**log_m_gr
    dr = 10**(log_dr)
    dm = 10**(log_dm)
    res = 1 / np.log(np.exp(1)+(a_gr/d)**n_gr)**m_gr
    res *= 1 - (np.log(1+dr/d)/np.log(1+dr/dm))**7
    return 100*res

def MSE(x, model, xdata, ydata):
    y_th = model(xdata, *x)
    residual = y_th - ydata
    MSE = 1/len(xdata) * np.sum(residual**2)
    return MSE

def quantile_loss(x, model, xdata, ydata, quantile):
    y_th = model(xdata, *x)
    residual = y_th - ydata
    L = np.where(residual < 0, - quantile * residual, - (quantile-1)* residual)
    return 1/len(xdata) * np.sum(L)

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
    bracket = [params[4],params[3]]
    res = scipy.optimize.root_scalar(lambda n: func(10**n, *params) - val, method="toms748", bracket=bracket)
    print(val, res)
    return float(10**res.root)


def fit(xdata, ydata):
    #from https://cdnsciencepub.com/doi/10.1139/t00-015, figure 11
    func = Fredlund_PSD
    bounds=[
        (-6,1), #log_a_gr_bounds
        (0,20), #log_n_gr_bounds
        (-1,5), #log_m_gr_bounds
        (-5,5), #log_dr_bounds
        (-6,0), #log_dm_bounds
    ]
    res = scipy.optimize.dual_annealing(MSE, bounds=bounds, args=(func, xdata, ydata), maxiter=1000, initial_temp=1e4, restart_temp_ratio=1e-3)
    return res

def fit_quantile(xdata, ydata, quantile, x0=None):
    func = Fredlund_PSD
    bounds=[
        (-6,1), #log_a_gr_bounds
        (0,20), #log_n_gr_bounds
        (-1,5), #log_m_gr_bounds
        (-5,5), #log_dr_bounds
        (-6,0), #log_dm_bounds
    ]
    #res = scipy.optimize.dual_annealing(quantile_loss, bounds=bounds, args=(func, xdata, ydata, quantile), maxiter=1000, initial_temp=1e4, restart_temp_ratio=1e-3)
    if x0 is None:
        res = scipy.optimize.differential_evolution(quantile_loss, bounds=bounds, args=(func, xdata, ydata, quantile), maxiter=1000, popsize=24)
    else:
        res = scipy.optimize.minimize(quantile_loss, x0=x0, bounds=bounds, args=(func, xdata, ydata, quantile))
    return res

def surface_specific(rho_s, func):
    # return specific surface based on the PSD
    Ss = 6/rho_s * scipy.integrate.quad()
    return Ss


###
# Prediction of saturated hydraulic conductivity
###  
# See also: https://link.springer.com/article/10.1007/s10064-012-0418-7

def predict_sat_perm_granular(porosity, Cu, D_10):
    # Predict the saturated permeability of the granular soil following Mbonimpa et al. (2002)
    # https://link.springer.com/article/10.1023/A:1016046214724
    # See eq. 13
    void_ratio = porosity / (1-porosity)
    gamma_w = 9807 #specific weight of water at 4° from https://www.engineeringtoolbox.com/water-density-specific-weight-d_595.html
    mu_w = 0.0015705 #water dynamic visocity at 4° from https://www.engineeringtoolbox.com/water-dynamic-kinematic-viscosity-d_596.html
    x = 2 #Mbonimpa et al. (2002)
    # K_g [same unit in input]^2
    #Note the formula coefficent X=1 m/cm, so X=100cm/cm
    K_g = 0.1 * 0.1 * void_ratio**(3+x) / (1+void_ratio) * Cu**(1/3) * D_10 * D_10
    return K_g


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
