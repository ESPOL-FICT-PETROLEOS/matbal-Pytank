from scipy.optimize import fsolve
# from pytank.aquifer import aquifer_carter_tracy,aquifer_fetkovich
from pytank.pvt_correlations import Bo_bw,comp_bw_nogas
from pytank.pvt_interp import interp_pvt_matbal
import numpy as np
from scipy.optimize import minimize
import matplotlib.pyplot as plt
from scipy import stats
import cvxpy as cp
import math
import pandas as pd
df_prod = pd.read_csv("p2_csv")
df_pvt= pd.read_csv("asd.csv")
ppvt_col = "P"
oil_fvf_col = "Bo"

def odia(n):
    return n * 365

def EBM(press, np, bo, cf, cw, sw0, boi, name, we,date,F,coeff):
    pi = 3215

    # Efw = boi * ((cf + cw * sw0) / (1 - sw0)) * (pi - press)
    # Eo = bo - boi
    # Et = Eo + Efw
    # # co = (bo - boi) / (boi * (pi - press))
    # # coeff = ((co * (1 - sw0)) + (cw * sw0) + cf) / (1 - sw0)
    # E= coeff * boi*(pi - press)
    # x = we / E
    # y = F / E
    y2 = press
    # Grafica
    # slope, intercept, r, p, se = stats.linregress(x, y)
    # N1 = slope
    # yt = intercept + (x * N1)
    # fig, ax = plt.subplots(figsize=(8, 6))
    # ax.scatter(x, y)
    # plt.plot(x, y)
    # ax.plot(x, yt, c='g', label='Regresion lineal')
    # print("Intercept(N)MMStb: \n", (intercept / 1e+6))
    # # Nombres de los ejes
    # text = " N [MMStb]: %.3f " % (intercept / 1000000)
    # plt.title('Gráfico Havlena y Odeh' + name)
    # plt.xlabel("WeBw/Eo+Efw")
    # plt.ylabel('F/Eo+Efw ')
    # plt.text(1e+7, 2e+8, text)
    # plt.legend()
    # plt.show()

    abajo = 654e+6 * boi * coeff
    abajo2 = 7.24078048e+08 * boi * coeff
    arriba = F - we
    #pr = pi - (arriba / abajo)
    #print(pr)
    pr2 = pi - (arriba / abajo2)
    # Grafica2
    fig, ax = plt.subplots(figsize=(20, 10))
    ax.scatter(np, y2, label='Presion Observada')
    #plt.plot(np, pr, c='g', label='Presion Calculada')
    plt.plot(np, date, c='r', label='Presion Optimizada')
    plt.title('Gráfico P vs t ', fontsize=15)
    plt.xlabel("Tiempo (Años)", fontsize=15)
    plt.ylabel('Presion (psia)', fontsize=15)
    #ax.set_ylim([0, 4000])
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    plt.legend(fontsize=15)
    ax.grid(axis='both', color='gray', linestyle='dashed')

    plt.show()

def mbal(p,pi,Np,wp,bo, cf, cw, sw0, boi,N,we,bw,wi):
    Eo = (bo - boi)
    Efw =  boi*(((cw * sw0) + cf) / (1 - sw0))*(pi-p)
    F = (Np * bo) + wp * bw
    return (N *(Eo+Efw) ) + (we * bw) + (wi*bw) - F
def aquifer_fetkovich(aq_radius, res_radius, aq_thickness, aq_por, ct, p, theta, k, water_visc,p_anterior,cum,pi):
    delta_t = 365
    wi = (math.pi / 5.615) * (aq_radius ** 2 - res_radius ** 2) * aq_thickness * aq_por
    f = theta / 360
    wei = ct * wi * pi * f
    rd = aq_radius / res_radius
    j = (0.00708 * k * aq_thickness * f) / (water_visc * (math.log(abs(rd))) )
    pa = pi*(1-(cum/wei))
    pr_avg = (p_anterior + p) / 2
    we = (wei / pi) * \
             (1 - np.exp((-1 * j * pi * delta_t) / wei)) * \
             (pa - pr_avg)
    cum_water_influx = cum + we
    return cum_water_influx

def press(p,Np,wp,cf,t,salinity,df_pvt,aq_radius,res_radius,aq_thickness, aq_por, theta, k, water_visc,p_anterior,cum,pi,sw0,N,wi):
    bo=interp_pvt_matbal(df_pvt,ppvt_col,oil_fvf_col,p)
    bw=Bo_bw(p,t,salinity,unit=1)
    boi=interp_pvt_matbal(df_pvt,ppvt_col,oil_fvf_col,pi)
    cw=comp_bw_nogas(p,t,salinity,unit=1)
    ct=cw+cf
    we = aquifer_fetkovich(aq_radius, res_radius, aq_thickness, aq_por, ct, p, theta, k, water_visc,p_anterior,cum,pi )
    return mbal(p,pi,Np,wp,bo, cf, cw, sw0, boi,N,we,bw,wi)

cf = 4.1e-6
t = 205
salinity= 55000
aq_radius = 35000
res_radius = 7000
aq_thickness= 100
aq_por=0.25
theta = 350
k = 43
cum=0
pi=3215
sw0 = 0.15
N= 654e+6
x0=3600
water_visc= np.exp(1.003-1.479e-2*(205)+1.982e-5*(205**2))

P_calculada = pd.DataFrame({'P_calculada': pi},index=[0])
for i in range(len(df_prod['P'])):
    Np=df_prod['Np'][i]
    wp = df_prod['Wp'][i]
    wi=df_prod['Wi'][i]
    p_anterior = P_calculada['P_calculada'][i]
    presion = fsolve(press,x0,args=(Np,wp,cf,t,salinity,df_pvt,aq_radius,res_radius,aq_thickness, aq_por, theta, k, water_visc,p_anterior,cum,pi,sw0,N,wi))
    x0 = presion[0]
    P_calculada = P_calculada._append({'P_calculada': presion[0]}, ignore_index=True)
    p_nueva = P_calculada['P_calculada'][i + 1]
    cw= comp_bw_nogas(p_nueva,t,salinity,unit=1)
    ct= cf + cw
    cum=aquifer_fetkovich(aq_radius, res_radius, aq_thickness, aq_por, ct, p_nueva, theta, k, water_visc,p_anterior,cum,pi)


def op(parametros,df,cf,t,salinity,df_pvt,sw0,cum,x0,k,water_visc,pi):
    N,aq_radius,theta=parametros
    P3 = pd.DataFrame({'P_calculada': pi}, index=[0])
    for i in range(len(df['P'])):
        Np = df_prod['Np'][i]
        wp = df_prod['Wp'][i]
        wi = df_prod['Wi'][i]
        p_anterior = P3['P_calculada'][i]
        presion = fsolve(press, x0, args=(
        Np, wp, cf, t, salinity, df_pvt, aq_radius, res_radius, aq_thickness, aq_por, theta, k, water_visc, p_anterior,
        cum, pi, sw0, N,wi))
        x0 = presion[0]
        P3 = P3._append({'P_calculada': presion[0]}, ignore_index=True)
        p_nueva = P3['P_calculada'][i + 1]
        cw = comp_bw_nogas(p_nueva, t, salinity, unit=1)
        ct = cf + cw
        cum = aquifer_fetkovich(aq_radius, res_radius, aq_thickness, aq_por, ct, p_nueva, theta, k, water_visc,
                                p_anterior, cum, pi)
    nueva_fila = pd.DataFrame({'Tank': 'tank_center', 'P': pi},
                              index=[0])
    df = pd.concat([nueva_fila, df]).reset_index(drop=True)
    df.iloc[0] = df.iloc[0].fillna(0)
    n_total=len(df['P'])
    error = (((df['P']-P3['P_calculada'])**2).sum())/n_total
    return error


df=df_prod

bnds=[(600e+6,700e+6),(0,None),(20,360)]
result = minimize(op, x0=[N,aq_radius,theta],args=(df,cf,t,salinity,df_pvt,sw0,cum,x0,k,water_visc,pi),bounds=bnds)
Nop,aq_radiusop, thetaop = result.x

print("Valor de sTOIIP que minimiza la función:", Nop)
print("Valor de aq radius que minimiza la función:", aq_radiusop/res_radius)
# print("Valor de res radiusque minimiza la función:", res_radiusop)
# print("Valor de aq_thickness que minimiza la función:", thicknessop)
print("Valor de theta que minimiza la función:",  thetaop)

N=Nop
aq_radius=aq_radiusop
#res_radius=res_radiusop
theta=thetaop
#aq_thickness=thicknessop

P4 = pd.DataFrame({'P_calculada': pi},index=[0])
for i in range(len(df_prod['P'])):
    Np = df_prod['Np'][i]
    wp = df_prod['Wp'][i]
    wi = df_prod['Wi'][i]
    p_anterior = P4['P_calculada'][i]
    presion = fsolve(press, x0, args=(
        Np, wp, cf, t, salinity, df_pvt, aq_radius, res_radius, aq_thickness, aq_por, theta, k, water_visc, p_anterior,
        cum, pi, sw0, N,wi))
    x0 = presion[0]
    P4 = P4._append({'P_calculada': presion[0]}, ignore_index=True)
    p_nueva = P4['P_calculada'][i + 1]
    cw = comp_bw_nogas(p_nueva, t, salinity, unit=1)
    ct = cf + cw
    cum = aquifer_fetkovich(aq_radius, res_radius, aq_thickness, aq_por, ct, p_nueva, theta, k, water_visc,
                            p_anterior, cum, pi)

nueva_fila = pd.DataFrame({'Date': '1987-09-01', 'P': pi, 'oil_fvf': 1.1},
                          index=[0])
df_prod = pd.concat([nueva_fila, df_prod]).reset_index(drop=True)
df_prod.iloc[0] = df_prod.iloc[0].fillna(0)
fig, ax = plt.subplots(figsize=(20, 10))
ax.scatter(df_prod['Np'], df_prod['P'], label='Presion Observada')
plt.plot(df_prod['Np'], P_calculada, c='g', label='Presion Calculada')
plt.plot(df_prod['Np'], P4, c='r', label='Presion Calculada')
plt.legend(fontsize=15)
plt.show()




