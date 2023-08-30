from pytank.aquifer import aquifer_carter_tracy
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats
from scipy.optimize import fsolve
from scipy.optimize import minimize
from scipy.optimize import least_squares
import cvxpy as cp
import math
from pytank.pvt_interp import interp_pvt_matbal
from pytank.pvt_correlations import Bo_bw, comp_bw_nogas

#%%

df_ta = pd.read_csv("mbal_Dataframe.csv")
df_ta2 = df_ta[df_ta["Tank"] == "tank_center"]
df_ta2.rename(
    columns={"PRESSURE_DATUM":"Pressure", "DATE": "Date"},inplace=True)
nueva_fila = pd.DataFrame({'Date': '1987-09-01', 'Pressure': 4500.00, 'oil_fvf': 1.1},
                          index=[0])
df_ta2 = pd.concat([nueva_fila, df_ta2]).reset_index(drop=True)
df_ta2.iloc[0] = df_ta2.iloc[0].fillna(0)
df_ta2['Date']=pd.to_datetime(df_ta2['Date'])
fecha= pd.to_datetime("1987-09-01")
df_ta2['time_setp'] = (df_ta2['Date']-fecha).dt.days

df_pvt = pd.read_csv("../tests/data_for_tests/full_example_1/pvt.csv")
df_pvt = df_pvt.fillna(method="ffill")
# time_step = time
ppvt_col = "Pressure"
oil_fvf_col = "Bo"
gas_fvf_col = "Bg"
gas_oil_rs_col = "GOR"
df_ta2['gas_prod_cum']= df_ta2['gas_prod_cum']*89
def mbal(p, pi, Np, wp, bo, cf, cw, sw0, boi, N, we,bw):
    #if rsi == rs:
    Eo = (bo - boi)
    Efw = boi * (((cw * sw0) + cf) / (1 - sw0)) * (pi - p)
    F = (Np * bo) + (wp * bw)
    funcion_P = (N * (Eo + Efw)) + (we * bw) - F
    # else:
    # Eo = (bo - boi) + (rsi-rs)*bg
    # Efw = (1+m) * boi * (((cw * sw0) + cf) / (1 - sw0)) * (pi - p)
    # Eg = boi*((bg/bgi)-1)
    # Gp= (gp-(Np*rs))*bg
    # F = (Np * bo) + (wp * bw) + Gp
    # funcion_P = (N * (Eo + Efw + m*Eg)) + (we * bw) - F

    return funcion_P


def carter(aq_por, ct, res_radius, aq_thickness, theta, aq_perm,
                         water_visc, pr, time,time_anterior,we,pi):
    pr_array = pr

    # Calculate the van Everdingen-Hurst water influx constant
    f = theta / 360
    b = 1.119 * aq_por * ct * (res_radius ** 2) * aq_thickness * f

    # Estimate dimensionless time (tD)
    cte = 0.006328 * aq_perm / (aq_por * water_visc * ct * (res_radius ** 2))
    td =  time * cte
    td2 = time_anterior*cte
    # Calculate the total pressure drop (Pi-Pn) as an array, for each time step n.
    pr_drop =  pi - pr_array

    # Estimate the dimensionless pressure
    pr_d =  0.5 * (np.log(td) + 0.80907)
    # Estimate the dimensionless pressure derivative
    e = 716.441 + (46.7984 * (td ** 0.5)) + (270.038 * td) + (71.0098 * (td ** 1.5))
    d = (1296.86 * (td ** 0.5)) + (1204.73 * td) + \
        (618.618 * (td ** 1.5)) + (538.072 * (td ** 2)) + \
        (142.41 * (td ** 2.5))
    pr_deriv = 1 / (2 * td)

    a1 = td - td2
    a2 = b * pr_drop
    a3 = we * pr_deriv
    a4 = pr_d
    a5 = td2 * pr_deriv
    cum_influx_water = we + (a1 * ((a2 - a3) / (a4 - a5)))
    we = cum_influx_water
    return we

def press(p, Np, wp, cf, t, salinity, df_pvt, res_radius, aq_thickness, aq_por, theta, k, water_visc,
          time,time_anterior, we, pi, sw0, N, boi):

    bo = interp_pvt_matbal(df_pvt, ppvt_col, oil_fvf_col, p)

    bw = Bo_bw(p, t, salinity, unit=1)
    cw = comp_bw_nogas(p, t, salinity, unit=1)
    ct = cw + cf
    we = carter(aq_por, ct, res_radius, aq_thickness, theta, k,
                         water_visc, p, time,time_anterior,we,pi)
    return mbal(p, pi, Np, wp, bo, cf, cw, sw0, boi, N, we,bw)


cf = 4.5e-6
t = 200
salinity = 30000
aq_radius = 15000
res_radius = 60
aq_thickness = 7
aq_por = 0.24
theta = 135
k = 42
water_visc = 0.35
cum = 0
pi = 4500
sw0 = 0.25
boi = interp_pvt_matbal(df_pvt, ppvt_col, oil_fvf_col, pi)
print(boi)
N = 70e+6
x0 = 3600
P_calculada = [pi]
x=len(df_ta2['Pressure'])


for i in range(x):
    if i !=33:
        Np = df_ta2['oil_prod_cum'][i+1]
        wp = df_ta2['water_prod_cum'][i+1]
        time = df_ta2['time_setp'][i+1]
        time_anterior= df_ta2['time_setp'][i]
        # Calculate current reservoir pressure given all other material balance variables through numeric solving.
        presion = fsolve(
            press, x0, args=(
                Np, wp, cf, t, salinity, df_pvt, res_radius, aq_thickness, aq_por, theta, k, water_visc,
                time, time_anterior, cum, pi, sw0, N, boi
            )
        )[0]
        print(f"Calculated Reservoir Pressure: {presion}")
        x0 = presion
        P_calculada.append(presion)
        cw = comp_bw_nogas(presion, t, salinity, unit=1)
        ct = cf + cw
        cum = carter(aq_por, ct, res_radius, aq_thickness, theta, k,
                             water_visc, x0, time,time_anterior,cum,pi)
        print(f"Calculated Cum: {cum}")
    else:
        print("termino")
    #print(f"Wei:{cum}")

# def op(parametros,df,cf,t,salinity,df_pvt,sw0,cum,x0,k,water_visc,pi):
#     N,aq_radius,res_radius,aq_thickness,theta=parametros
#     P3 = pd.DataFrame({'P_calculada': pi}, index=[0])
#     for i in range(len(df['DATE'])):
#         Np = df['oil_prod_cum'][i]
#         wp = df['water_prod_cum'][i]
#         p_anterior = P3['P_calculada'][i]
#         presion = fsolve(press, x0, args=(
#         Np, wp, cf, t, salinity, df_pvt, aq_radius, res_radius, aq_thickness, aq_por, theta, k, water_visc, p_anterior,
#         cum, pi, sw0, N))
#         x0 = presion[0]
#         P3 = P3._append({'P_calculada': presion[0]}, ignore_index=True)
#         p_nueva = P3['P_calculada'][i + 1]
#         cw = comp_bw_nogas(p_nueva, t, salinity, unit=1)
#         ct = cf + cw
#         cum = aquifer_fetkovich(aq_radius, res_radius, aq_thickness, aq_por, ct, p_nueva, theta, k, water_visc,
#                                 p_anterior, cum, pi)
#     nueva_fila = pd.DataFrame({'Tank': 'tank_center', 'Pressure': 3700.00},
#                               index=[0])
#     df = pd.concat([nueva_fila, df]).reset_index(drop=True)
#     df.iloc[0] = df.iloc[0].fillna(0)
#     n_total=len(df['Pressure'])
#     error = (((df['Pressure']-P3['P_calculada'])**2).sum())/n_total
#     return error
#
#
# df=df_ta2
#
# bnds=[(1e+6,100e+6),(0,None),(0,None),(0,None),(20,360)]
# result = minimize(op, x0=[N,aq_radius,res_radius,aq_thickness,theta],args=(df,cf,t,salinity,df_pvt,sw0,cum,x0,k,water_visc,pi),bounds=bnds)
# Nop,aq_radiusop,res_radiusop,thicknessop, thetaop = result.x
#
# print("Valor de sTOIIP que minimiza la función:", Nop)
# print("Valor de aq radius que minimiza la función:", aq_radiusop)
# print("Valor de res radiusque minimiza la función:", res_radiusop)
# print("Valor de aq_thickness que minimiza la función:", thicknessop)
# print("Valor de theta que minimiza la función:",  thetaop)
#
# N=Nop
# aq_radius=aq_radiusop
# res_radius=res_radiusop
# theta=thetaop
# aq_thickness=thicknessop
#
# P4 = pd.DataFrame({'P_calculada': 3700.00},index=[0])
# for i in range(len(df_ta2['DATE'])):
#     Np = df_ta2['oil_prod_cum'][i]
#     wp = df_ta2['water_prod_cum'][i]
#     p_anterior = P4['P_calculada'][i]
#     presion = fsolve(press, x0, args=(
#         Np, wp, cf, t, salinity, df_pvt, aq_radius, res_radius, aq_thickness, aq_por, theta, k, water_visc, p_anterior,
#         cum, pi, sw0, N))
#     x0 = presion[0]
#     P4 = P4._append({'P_calculada': presion[0]}, ignore_index=True)
#     p_nueva = P4['P_calculada'][i + 1]
#     cw = comp_bw_nogas(p_nueva, t, salinity, unit=1)
#     ct = cf + cw
#     cum = aquifer_fetkovich(aq_radius, res_radius, aq_thickness, aq_por, ct, p_nueva, theta, k, water_visc,
#                             p_anterior, cum, pi)
#%%

df_ta2["Date"] = pd.to_datetime(df_ta2["Date"])
df_ta2.iloc[0] = df_ta2.iloc[0].fillna(0)
#%%
fig, ax = plt.subplots(figsize=(15, 10))
ax.scatter(df_ta2['Date'].dt.year, df_ta2['Pressure'], label='Presion Observada')
plt.plot(df_ta2['Date'].dt.year, P_calculada, c='g', label='Presion Calculada')
# plt.plot(df_ta2['Date'], P4, c='r', label='Presion Calculada')
plt.title('Gráfico P vs t ', fontsize=15)
plt.xlabel("Tiempo (Años)", fontsize=15)
plt.ylabel('Presion (psia)', fontsize=15)
#ax.set_ylim([400, 4000])
# plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
ax.grid(axis='both', color='gray', linestyle='dashed')
plt.legend(fontsize=15)
plt.show()




