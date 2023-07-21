from pytank.aquifer import aquifer_fetkovich,aquifer_carter_tracy
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats
from scipy.optimize import fsolve
from scipy.optimize import minimize
import math
from pytank.pvt_interp import interp_pvt_matbal
from pytank.pvt_correlations import Bo_bw, comp_bw_nogas

def odia(n):
    return n * 30
def F(np, bo, wp):
    F = (np * bo) + wp
    return F


def coeff(p, sw0, cw, cf, bo, boi):
    pi = 3700
    co = (bo - boi) / (boi * (pi - p))
    coeff = ((co * (1 - sw0)) + (cw * sw0) + cf) / (1 - sw0)
    return coeff
def EBM(press, np, wp, bo, cf, cw, sw0, boi, name, we,date):
    pi = 3700
    F = (np * bo) + wp
    #Efw = boi * ((cf + cw * sw0) / (1 - sw0)) * (pi - press)
    #Eo = bo - boi
    co = (bo - boi) / (boi * (pi - press))
    coeff = ((co * (1 - sw0)) + (cw * sw0) + cf) / (1 - sw0)
    Et = coeff*boi*(pi - press)
    x = np
    y = F / Et

    y2 = press
    # Grafica
    # slope, intercept, r, p, se = stats.linregress(x, y)
    # N1 = slope
    # yt = intercept + (x * N1)
    fig, ax = plt.subplots(figsize=(8, 6))
    ax.scatter(date, y)

    #ax.plot(x, yt, c='g', label='Regresion lineal')
    #print("Intercept(N)MMStb: \n", (intercept / 1e+6))

    # Nombres de los ejes
    #text = " N [MMStb]: %.3f " % (intercept / 1000000)
    plt.title('Gráfico Campbell')
    plt.xlabel("Tiempo (Años)")
    plt.ylabel('F/Eo+Efw ')
    #plt.text(1e+7, 2e+8, text)
    #plt.legend()
    plt.show()

    # # pr
    # abajo = intercept * boi * coeff
    # abajo2 = 19543283.851659186 * boi * coeff
    # arriba = F - we
    # pr = pi - (arriba / abajo)
    # pr2 = pi - (arriba / abajo2)
    #
    # # Grafica2
    # fig, ax = plt.subplots(figsize=(20, 10))
    # ax.scatter(date, y2, label='Presion Observada')
    # plt.plot(date, pr, c='g', label='Presion Calculada')
    # plt.plot(date, p2, c='r', label='Presion Optimizada')
    # plt.title('Gráfico P vs t ', fontsize=15)
    # plt.xlabel("Tiempo (Años)", fontsize=15)
    # plt.ylabel('Presion (psia)', fontsize=15)
    # ax.set_ylim([0,4000])
    # plt.xticks(fontsize=15)
    # plt.yticks(fontsize=15)
    # plt.legend(fontsize=15)
    # ax.grid(axis='both', color='gray', linestyle='dashed')
    #
    # plt.show()
df_pvt = pd.read_csv("../tests/data_for_tests/full_example_1/pvt.csv")
df_pvt = df_pvt.fillna(method="ffill")
ppvt_col = "Pressure"
oil_fvf_col = "Bo"
nueva_fila = pd.DataFrame({'Date': '1981-08-01', 'Tank': 'tank_center', 'Pressure': 3700.00},
                          index=[0])
df_ta = pd.read_csv("mbal_Dataframe3.csv")
df_center = df_ta[df_ta["Tank"] == "tank_center"]
#df_center = df_center.loc[(df_center['DATE'] >= '2002-09-01') & (df_center['DATE'] < '2007-09-01')]

df_center = pd.concat([nueva_fila, df_center]).reset_index(drop=True)
df_center.iloc[0] = df_center.iloc[0].fillna(0)
date = df_center['Date']
# %% Calcular Influjo de agua
time = []
for i in range(len(df_center['Date'])):
    time.append(odia(i))
p = df_center['Pressure']
aq_radius = 14000
res_radius = 2000
aq_thickness = 20
phi = 0.25
ct = 0.000008
pr = p.to_numpy()
theta = 300
k = 20
water_visc = 0.4
time_step = time
w = aquifer_fetkovich(aq_radius, res_radius, aq_thickness, phi, ct, pr, theta, k, water_visc, time_step,
                      boundary_type='no_flow', flow_type='radial', width=None, length=None)
# w=aquifer_carter_tracy(aq_por, ct, res_radius, aq_thickness, theta, aq_perm,
#                           water_visc, pr, time)
# %% Adjuntar el dataframe del influjo de agua al dataframe principal
#df_center['AWe'] = df_center.get("AWe", w['Delta We'].to_list())
df_center['Cwe'] = df_center.get("Cwe", w['Cumulative We'].to_list())
# D = 'START_DATETIME'
# df_center = df_center.drop([D], axis=1)
#
# %% Datos para el balance de materiales
df_center = df_center.drop(index=0)
press = df_center['Pressure']
np = df_center['oil_prod_cum']
wp = df_center['water_prod_cum']
gp = df_center['gas_prod_cum']
bo = df_center['oil_fvf']
bg = df_center['gas_fvf']
rs = df_center['gas_oil_rs_col']
we = df_center['Cwe']
name = "Tank center pressure avg"
df_center["Date"] = pd.to_datetime(df_center["Date"])
date=df_center['Date']
sw0 = 0.25
cf = 5e-6
cw = 3e-6
boi=1.1040
#boi = interp_pvt_matbal(df_pvt, ppvt_col, oil_fvf_col, 3700)
coeff = coeff(p, sw0, cw, cf, bo, boi)
F = F(np, bo, wp)
df_center['F'] = F
df_center['coeff'] = coeff
EBM(press, np, wp, bo, cf, cw, sw0, boi, name, we,date)
# def f(we,df):
#     return (1.75e+8 * 1.11 * df['coeff'] * (3700 - df['Pressure'])) + we - df['F']
#
# def funcion(P,df):
#     return df['F']-(1.75e+8*1.11*df['coeff']*(3700-P))-df['Cwe']
# P2 = pd.DataFrame()
# wei=pd.DataFrame()
# x0=0
# for i in range(len(date)):
#     fila=df_center.iloc[i]
#     df=fila
#     r=fsolve(f,x0,args=df)
#     print(r)
#     x0=r[0]
#     wei=wei._append( {'wei': r[0]},ignore_index=True)
# df_center['Cwe2'] = df_center.get("Cwe2", wei['wei'].to_list())
# x1=3700
# for i in range(len(date)):
#     fila=df_center.iloc[i]
#     df=fila
#     re=fsolve(funcion,x1,args=df)
#     x1=re[0]
#     P2=P2._append( {'P2': re[0]},ignore_index=True)
#
# print(df_center['Cwe'])
# print(df_center['Cwe2'])

# p2=P2
# p2=p2.drop(index=0)
# print(P2)






# def f(para,df):
#     delta_t = 365
#     pr_array = df['Pressure'].to_numpy()
#     aq_radius, res_radius,k,t = para
#     wi = (math.pi / 5.615) * (aq_radius ** 2 - res_radius ** 2) * 50 * 0.25
#     f = t / 360
#     wei = 8e-6 * wi * 3700 * f
#     rd = aq_radius / res_radius
#     j = (0.00708 * k * 50 * f) / (0.8 * (math.log(abs(rd)))-0.75)
#     cum_water_influx = 0
#     pr = 3700
#     pa = 3700
#     df2 = pd.DataFrame(columns=['Delta We'])
#     for ip in range(len(pr_array)):
#         pr_avg = (pr + pr_array[ip]) / 2
#         we = (wei / pr_array[0]) * \
#                 (1 - math.exp((-1 * j * pr_array[0] * delta_t) / wei)) * \
#                 (pa - pr_avg)
#
#         pr = pr_array[ip]
#         cum_water_influx = cum_water_influx + we
#         pa = pr_array[0] * (1 - (cum_water_influx / wei))
#         df2 = df2._append({'Delta We': we, 'Cumulative We': cum_water_influx}, ignore_index=True)
#     error=(df['Cwe2']-df2['Cumulative We'])**2
#     error2=math.sqrt((error.sum()/len(df['Date'])))
#     return error2
# x0 = [18000,2000,40,45]
# result = minimize(f, x0,args=df)
# x_min, y_min,z,t = result.x
# print("Valor de aq radius que minimiza la función:", x_min)
# print("Valor de res radius que minimiza la función:", y_min)
# print("Valor de k que minimiza la función:", z)
# print("Valor de Theta que minimiza la función:", t)
#
# def funcion_objetivo(parametros, df):
#     N= parametros
#     presion_calculada = 3700 - (
#                 (df['F'] - df['Cwe']) / (N * 1.11 * df['coeff']))
#     diferencia = df['Pressure'] - presion_calculada
#     df['d'] = diferencia ** 2
#     error_cuadratico = math.sqrt(df['d'].sum()/len(df['F']))
#     return error_cuadratico
# x1=[80e+6]
# r=minimize(funcion_objetivo, x1,args=df)
# n = r.x
# print("Valor de n que minimiza la función:", n)
# print("Valor de boi que minimiza la función:", boi_op)
#EBM(press, np, bo, cf, cw, sw0, boi, name, we,date,p2)