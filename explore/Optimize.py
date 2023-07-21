from pytank.aquifer import aquifer_fetkovich
import pandas as pd
import cvxpy as cp
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
    k = 20
    water_visc = 0.75
    pi = 3676
    aq_radius=9455
    res_radius=600
    time = []
    for i in range(len(df['DATE'])):
        time.append(odia(i))
    time_step = time
    pr = df['PRESSURE_DATUM'].to_numpy()
    N = parametros
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
df_center = df_ta[df_ta["Tank"] == "tank_center"]
df_center = df_center.loc[(df_center['DATE'] >= '2002-09-01') & (df_center['DATE'] < '2008-09-01')]
df_center = pd.concat([nueva_fila, df_center]).reset_index(drop=True)
df_center.iloc[0] = df_center.iloc[0].fillna(0)
date = df_center['DATE']
time = []
for i in range(len(df_center['DATE'])):
    time.append(odia(i))
time_step = time
# Datos
p = df_center['PRESSURE_DATUM']
np = df_center['oil_prod_cum']
wp = df_center['water_prod_cum']
gp = df_center['gas_prod_cum']
bo = df_center['oil_fvf']
bg = df_center['gas_fvf']
rs = df_center['gas_oil_rs_col']
date = df_center['DATE']
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
df_center['F'] = F
df_center['coeff'] = coeff
df_center['time_step'] = time_step
df= df_center
# %% Calcular el error cuadratico
# parametros = 18253409.514814563, 9455, 600
# df = df_center
# e = funcion_objetivo(parametros, df)
# print("error cuadratico:", e)
# %% Calcular los valores optimos

#valores_iniciales = 18253409.514814563, 9455, 600

# resultado = minimize(funcion_objetivo, valores_iniciales, args=df)
# print(resultado,"hola")
# # Obtención de los valores óptimos de los parámetros
# N,raquier_optimo, res_optimo = resultado.x
# print("Radio del acuífero óptimo:", resultado.x)

N = cp.Variable()


# Define the objective function
objective = cp.Minimize(funcion_objetivo([N], df))

# Define the constraints (if any)
# For example, if you want to add constraints on the variables, you can do it here.

# Define the optimization problem
problem = cp.Problem(objective)

# Solve the problem
result = problem.solve()

# Get the optimized values
optimized_N = N.value


print("Optimized Values:")
print("N:", optimized_N)
