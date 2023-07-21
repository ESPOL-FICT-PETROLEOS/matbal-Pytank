from pytank.aquifer import aquifer_fetkovich
import pandas as pd

from scipy.optimize import minimize


def odia(n):
    return n * 365

def F(np, bo, wp):
    F = (np * bo) + wp
    return F


def coeff(p, sw0, cw, cf, bo, boi):
    pi = 3676
    co = (bo - boi) / (boi * (pi - p))
    coeff = ((co * (1 - sw0)) + (cw * sw0) + cf) / (1 - sw0)
    return coeff


def funcion_objetivo(parametros, df):
    boi = 0.86
    aq_thickness = 50
    phi = 0.26
    ct = 0.000007
    theta = 120
    k = 60
    water_visc = 0.75
    pi = 3676
    pr = df['PRESSURE_DATUM'].to_numpy()
    N, aq_radius, res_radius = parametros
    w = aquifer_fetkovich(aq_radius, res_radius, aq_thickness, phi, ct, pr, theta, k, water_visc, time_step,
                          boundary_type='no_flow', flow_type='radial', width=None, length=None)
    df['CWe'] = df.get("CWe", w['Cumulative We'].to_list())
    df = df.drop(index=0)
    presion_calculada = pi - (
                (df['F'] - df['CWe']) / (N * boi * df['coeff']))  # Reemplaza Pr con tu función de balance de materia
    diferencia = df['PRESSURE_DATUM'] - presion_calculada
    df['d'] = diferencia ** 2
    error_cuadratico = df['d'].sum()
    return error_cuadratico

# %% Dataframe
nueva_fila = pd.DataFrame({'DATE': '1981-08-01', 'Tank': 'tank_center', 'PRESSURE_DATUM': 3676.00, 'oil_fvf': 1.1},
                          index=[0])
df_ta = pd.read_csv("mbal_Dataframe2.csv")
df_tanknorth2 = df_ta[df_ta["Tank"] == "tank_center"]
df_tanknorth2 = df_tanknorth2.loc[(df_tanknorth2['DATE'] >= '2002-09-01') & (df_tanknorth2['DATE'] < '2008-09-01')]
df_tanknorth2 = pd.concat([nueva_fila, df_tanknorth2]).reset_index(drop=True)
df_tanknorth2.iloc[0] = df_tanknorth2.iloc[0].fillna(0)
date = df_tanknorth2['DATE']
time = []
for i in range(len(df_tanknorth2['DATE'])):
    time.append(odia(i))
time_step = time
# Datos
p = df_tanknorth2['PRESSURE_DATUM']
np = df_tanknorth2['oil_prod_cum']
wp = df_tanknorth2['water_prod_cum']
gp = df_tanknorth2['gas_prod_cum']
bo = df_tanknorth2['oil_fvf']
bg = df_tanknorth2['gas_fvf']
rs = df_tanknorth2['gas_oil_rs_col']
date = df_tanknorth2['DATE']
name = "Tank center pressure avg"
sw0 = 0.15
cf = 3.5e-6
cw = 3.62e-6
boi = 0.86
# %% Calcular coeff (compresibilidad efectiva)
coeff = coeff(p, sw0, cw, cf, bo, boi)
# %% Calcular F
F = F(np, bo, wp)
# %% Ajuntar estos valores al dataframe
df_tanknorth2['F'] = F
df_tanknorth2['coeff'] = coeff
df_tanknorth2['time_step'] = time_step
# %% Calcular el error cuadratico
parametros = 16238924.341140533, 9455, 600
df = df_tanknorth2
e = funcion_objetivo(parametros, df)
print(e)
# %% Calcular los valores optimos

valores_iniciales = 16238924.341140533, 9455, 600

resultado = minimize(funcion_objetivo, valores_iniciales, args=df)

# Obtención de los valores óptimos de los parámetros
N_optimo, raquier_optimo, rreservorio_optimo = resultado.x
print(
    "N:" + str(N_optimo) + "\nr del acuifero:" + str(raquier_optimo) + "\nr del reservorio:" + str(rreservorio_optimo))
