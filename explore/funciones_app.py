from scipy import stats
import matplotlib.pyplot as plt
import pandas as pd
import math
from pytank.pvt_interp import interp_pvt_matbal
from pytank.pvt_correlations import Bo_bw, comp_bw_nogas
import numpy as np
from scipy.optimize import fsolve
def bw(presion,t,salinity):
    Bw=[]
    for i in range(len(presion)):
        p=presion[i]
        bw = Bo_bw(p, t, salinity, unit=1)
        Bw.append(bw)
    Bw = pd.DataFrame({"Bw": Bw})
    return Bw

def cw(presion,t,salinity):
    Cw=[]
    for i in range(len(presion)):
        p=presion[i]
        cw = comp_bw_nogas(p, t, salinity, unit=1)
        Cw.append(cw)
    Cw = pd.DataFrame({"Cw": Cw})
    return Cw

def Campbell(p, np,wp, bo, cf, sw0, boi, date,pi,t,salinity):
    Bw=bw(p, t, salinity)['Bw']
    Cw=cw(p, t, salinity)['Cw']
    Eo = (bo - boi)
    Efw = boi * (((Cw * sw0) + cf) / (1 - sw0)) * (pi - p)
    F = (np * bo) + (wp*Bw)
    y = F/(Eo+Efw)
    x = date
    data = pd.DataFrame({'Date': x, 'F/Eo+Efw': y})
    return data

def G_method(pr, np,wp, bo, cf, sw0, boi, we,pi,t,salinity):
    Bw = bw(pr, t, salinity)['Bw']
    Cw = cw(pr, t, salinity)['Cw']
    Eo = (bo - boi)
    Efw = boi * (((Cw * sw0) + cf) / (1 - sw0)) * (pi - pr)
    F = (np * bo) + (wp * Bw)
    y = F / (Eo + Efw)
    x = (we*Bw)/ (Eo + Efw)
    data = pd.DataFrame({'We*Bw/Et': x, 'F/Eo+Efw': y})
    return data



def EBM(p, pi, Np, wp, bo, cf, cw, sw0, boi, N, we,bw):
    Eo = (bo - boi)
    Efw = boi * (((cw * sw0) + cf) / (1 - sw0)) * (pi - p)
    F = (Np * bo) + (wp * bw)
    funcion_P = (N * (Eo + Efw)) + (we * bw) - F
    return funcion_P

def aquifer_fetkovich(aq_radius, res_radius, aq_thickness, aq_por, ct, p, theta, k, water_visc, p_anterior, cum, pi):
    delta_t = 365
    wi = (math.pi / 5.615) * (aq_radius ** 2 - res_radius ** 2) * aq_thickness * aq_por
    f = theta / 360
    wei = ct * wi * pi * f
    rd = aq_radius / res_radius
    j = (0.00708 * k * aq_thickness * f) / (water_visc * (math.log(abs(rd))))
    pa = pi * (1 - (cum / wei))
    pr_avg = (p_anterior + p) / 2
    we = (wei / pi) * \
         (1 - np.exp((-1 * j * pi * delta_t) / wei)) * \
         (pa - pr_avg)
    cum_water_influx = cum + we
    return cum_water_influx

def press(p, Np, wp, cf, t, salinity, df_pvt, aq_radius, res_radius, aq_thickness, aq_por, theta, k, water_visc,
          p_anterior, cum, pi, sw0, N, boi, ppvt_col, oil_fvf_col):
    #Parametros que depende de la presión
    bo = interp_pvt_matbal(df_pvt, ppvt_col, oil_fvf_col, p)
    bw = Bo_bw(p, t, salinity, unit=1)
    cw = comp_bw_nogas(p, t, salinity, unit=1)
    ct = cw + cf
    we = aquifer_fetkovich(aq_radius, res_radius, aq_thickness, aq_por, ct, p, theta, k, water_visc, p_anterior, cum,
                           pi)
    return EBM(p, pi, Np, wp, bo, cf, cw, sw0, boi, N, we,bw)

def calcuted_pressure(Np_frame, wp_frame, cf, t, salinity, df_pvt, aq_radius, res_radius, aq_thickness, aq_por, theta, k, water_visc,
            pi, sw0, N, ppvt_col, oil_fvf_col):
    #Valores iniciales
    boi = interp_pvt_matbal(df_pvt, ppvt_col, oil_fvf_col, pi)
    cum = 0
    P_calculada = [pi]
    x0=pi
    # Iteración de cada uno de los años
    for i in range(len(Np_frame)):
        Np = Np_frame[i]
        wp = wp_frame[i]
        p_anterior = P_calculada[i]
        # Calculate current reservoir pressure given all other material balance variables through numeric solving.
        presion = fsolve(
            press, x0, args=(
                Np, wp, cf, t, salinity, df_pvt, aq_radius, res_radius, aq_thickness, aq_por, theta, k, water_visc,
                p_anterior, cum, pi, sw0, N, boi, ppvt_col, oil_fvf_col
            )
        )[0]
        x0 = presion
        P_calculada.append(presion)
        cw = comp_bw_nogas(presion, t, salinity, unit=1)
        ct = cf + cw
        cum = aquifer_fetkovich(aq_radius, res_radius, aq_thickness, aq_por, ct, presion, theta, k, water_visc, p_anterior,
                                cum, pi)
    return P_calculada