from pytank.aquifer import aquifer_fetkovich
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats


def odia(n):
    return n * 365

def EBM(press, np, wp, bo, cf, cw, sw0, boi, name, we, date):
    pi = 3676
    F = (np * bo) + wp
    Efw = boi * ((cf + cw * sw0) / (1 - sw0)) * (pi - press)
    Eo = bo - boi
    Et = Eo + Efw
    x = we / Et
    y = F / Et
    co = (bo - boi) / (boi * (pi - press))
    coeff = ((co * (1 - sw0)) + (cw * sw0) + cf) / (1 - sw0)
    y2 = press
    # Grafica
    slope, intercept, r, p, se = stats.linregress(x, y)
    N1 = slope
    yt = intercept + (x * N1)
    fig, ax = plt.subplots(figsize=(8, 6))
    ax.scatter(x, y)
    plt.plot(x, y)
    ax.plot(x, yt, c='g', label='Regresion lineal')
    print("Intercept(N)MMStb: \n", (intercept))

    # Nombres de los ejes
    text = " N [MMStb]: %.3f " % (intercept / 1000000)
    plt.title('Gráfico Havlena y Odeh' + name)
    plt.xlabel("WeBw/Eo+Efw")
    plt.ylabel('F/Eo+Efw ')
    plt.text(1e+7, 2e+7, text)
    plt.legend()
    plt.show()

    # pr
    abajo = intercept * boi * coeff
    abajo2 = 15956762.922427185 * boi * coeff
    arriba = F - we
    pr = pi - (arriba / abajo)
    pr2 = pi - (arriba / abajo2)

    # Grafica2
    fig, ax = plt.subplots(figsize=(15, 10))
    ax.scatter(date, y2, label='Presion Observada')
    plt.plot(date, pr, c='g', label='Presion Calculada')
    plt.plot(date, pr2, c='r', label='Presion Optimizada')
    plt.title('Gráfico P vs t ', fontsize=15)
    plt.xlabel("Tiempo (Años)", fontsize=15)
    plt.ylabel('Presion (psia)', fontsize=15)
    ax.set_ylim([600,2400])
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    plt.legend(fontsize=15)
    ax.grid(axis='both', color='gray', linestyle='dashed')

    plt.show()


nueva_fila = pd.DataFrame({'DATE': '1981-08-01', 'Tank': 'tank_center', 'PRESSURE_DATUM': 3676.00, 'oil_fvf': 1.1},
                          index=[0])
df_ta = pd.read_csv("mbal_Dataframe2.csv")
df_tanknorth2 = df_ta[df_ta["Tank"] == "tank_center"]
df_tanknorth2 = df_tanknorth2.loc[(df_tanknorth2['DATE'] >= '2002-09-01') & (df_tanknorth2['DATE'] < '2008-09-01')]
df_tanknorth2 = pd.concat([nueva_fila, df_tanknorth2]).reset_index(drop=True)
df_tanknorth2.iloc[0] = df_tanknorth2.iloc[0].fillna(0)
date = df_tanknorth2['DATE']
# %% Calcular Influjo de agua
time = []
for i in range(len(df_tanknorth2['DATE'])):
    time.append(odia(i))
p = df_tanknorth2['PRESSURE_DATUM']
aq_radius = 9455
res_radius = 600
aq_thickness = 50
phi = 0.26
ct = 0.000007
pr = p.to_numpy()
theta = 120
k = 60
water_visc = 0.75
time_step = time
w = aquifer_fetkovich(aq_radius, res_radius, aq_thickness, phi, ct, pr, theta, k, water_visc, time_step,
                      boundary_type='no_flow', flow_type='radial', width=None, length=None)
# %% Adjuntar el dataframe del influjo de agua al dataframe principal
df_tanknorth2['AWe'] = df_tanknorth2.get("AWe", w['Delta We'].to_list())
df_tanknorth2['Cwe'] = df_tanknorth2.get("Cwe", w['Cumulative We'].to_list())
D = 'START_DATETIME'
df_tanknorth2 = df_tanknorth2.drop([D], axis=1)
df_tanknorth2 = df_tanknorth2.drop(index=0)
# %% Datos para el balance de materiales
press = df_tanknorth2['PRESSURE_DATUM']
np = df_tanknorth2['oil_prod_cum']
wp = df_tanknorth2['water_prod_cum']
gp = df_tanknorth2['gas_prod_cum']
bo = df_tanknorth2['oil_fvf']
bg = df_tanknorth2['gas_fvf']
rs = df_tanknorth2['gas_oil_rs_col']
date = df_tanknorth2['DATE']
we = df_tanknorth2['Cwe']
name = "Tank center pressure avg"
sw0 = 0.15
cf = 3.5e-6
cw = 3.62e-6
boi = 0.86
EBM(press, np, wp, bo, cf, cw, sw0, boi, name, we, date)
