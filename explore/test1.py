#%%
from pytank.material_balance import campbell_function, ho_terms_equation, \
    havlena_and_odeh
from pytank.aquifer import aquifer_fetkovich
import pandas as pd
import matplotlib.pyplot as plt
def odia(n):
    return(n*365)
from scipy import stats
def EBM(p,np,wp,bo,cf,cw,sw0,boi, name,date):

    pi = 3676
    F = (np*bo)+wp
    Efw = boi*((cf + cw * sw0) / (1 - sw0)) * (pi - p)
    Eo = bo-boi
    Et = Eo+Efw
    x = date
    y = F / Et
    y2 = p
    print(p,F,Eo)
    #Grafica
    #slope, intercept, r, p, se = stats.linregress(x, y) #encontrar la pendiente
    #N1=slope
    #yt = intercept + (x * N1)
    fig, ax = plt.subplots(figsize=(8, 6))
    ax.scatter(x, y)
    plt.plot(x, y)
    #ax.plot(x, yt, c='g',label='Hopecreek')
    #print("pendiente(N): \n",intercept)

    #Nombres de los ejes
    #text = " N [MMStb]: %.3f " % intercept
    plt.title('Gráfico Havlena y Odeh' + name)
    plt.xlabel("WeBw/Eo+Efw")
    plt.ylabel('F/Eo+Efw ')
    #plt.text(0.28, 1750000, text)
    plt.show()

    # pr
    pr = pi - (((F - intercept * Eo) - we)) / ((cf + cw * sw0 / (1 - sw0)) * boi * intercept)
    print(pr)
    # Grafica2
    fig, ax = plt.subplots(figsize=(15, 10))
    ax.scatter(date, y2)
    plt.plot(date, y2)
    ax.scatter(date, pr)
    plt.plot(date, pr)

    plt.show()

#%% Load the mbal dataframe
df_tank = pd.read_csv("mbal_Dataframe.csv")
df_ta=pd.read_csv("mbal_Dataframe2.csv")
D='START_DATETIME'
df_ta=df_ta.drop([D], axis=1)
#%% Creation of a mbal dataframe for each tank
df_tankcenter = df_tank[df_tank["Tank"] == "tank_center"]
df_tankcenter2 = df_ta[df_ta["Tank"] == "tank_center"]
df_tanknorth2 = df_ta[df_ta["Tank"] == "tank_north"]
#%% Additional numerical data
sw0 = 0.15
cf = 3.5e-6
cw = 3.62e-6
boi = 0.86
aq_radius = 46000
res_radius = 9200
aq_thickness = 100
phi = 0.25
ct = 0.000007
theta = 140
k = 200
water_visc = 0.55

menu=int(input("Menu:h  \n 1.tank_north(avg) \n 2. tank_center \n 3.tank_center(avg)"))
while menu !=0:
    if menu ==1:
        name = "Tank north 1989-2003"
        #df_tanknorth2 = df_tanknorth2.loc[(df_tanknorth2['DATE'] >= '1989-08-01') & (df_tanknorth2['DATE'] < '2003-08-01')]
        date = df_tanknorth2['DATE']
        time = []
        for i in range(len(df_tanknorth2['DATE']) ):
            time.append(odia(i))
        p = df_tanknorth2['PRESSURE_DATUM']
        np = df_tanknorth2['oil_prod_cum']
        wp = df_tanknorth2['water_prod_cum']
        gp = df_tanknorth2['gas_prod_cum']
        bo = df_tanknorth2['oil_fvf']
        bg = df_tanknorth2['gas_fvf']
        rs = df_tanknorth2['gas_oil_rs_col']
        pr = p.to_numpy()
        time_step = time
        #w=aquifer_fetkovich(aq_radius, res_radius, aq_thickness, phi, ct, pr,theta, k, water_visc, time_step,boundary_type='no_flow', flow_type='radial',width=None, length=None)
        #df_tanknorth2['AWe'] = df_tanknorth2.get("AWe", w['Delta We'].to_list())
        #df_tanknorth2['Cwe'] = df_tanknorth2.get("Cwe", w['Cumulative We'].to_list())
        #we = df_tanknorth2['Cwe']
        EBM(p,np,wp,bo,cf,cw,sw0,boi, name,date)

    elif menu==2:
        name = "Tank center"
        p = df_tankcenter['PRESSURE_DATUM']
        np = df_tankcenter['oil_prod_cum']
        wp = df_tankcenter['water_prod_cum']
        gp = df_tankcenter['gas_prod_cum']
        bo = df_tankcenter['oil_fvf']
        bg = df_tankcenter['gas_fvf']
        rs = df_tankcenter['gas_oil_rs_col']
        EBM(p,np,wp,bo,bg,cf,cw,sw0,rs,boi,name,date,we)
    elif menu==3:
        name = "Tank center 1988-1995 "
        #df_tankcenter2 = df_tankcenter2.loc[
            #(df_tankcenter2['DATE'] >= '2002-09-01') & (df_tankcenter2['DATE'] < '2008-09-01')]
        date = df_tankcenter2['DATE']
        time = []
        for i in range(len(df_tankcenter2['DATE'])):
            time.append(odia(i))
        p = df_tankcenter2['PRESSURE_DATUM']
        np = df_tankcenter2['oil_prod_cum']
        wp = df_tankcenter2['water_prod_cum']
        gp = df_tankcenter2['gas_prod_cum']
        bo = df_tankcenter2['oil_fvf']
        bg = df_tankcenter2['gas_fvf']
        rs = df_tankcenter2['gas_oil_rs_col']
        pr = p.to_numpy()
        time_step = time
        w = aquifer_fetkovich(aq_radius, res_radius, aq_thickness, phi, ct, pr, theta, k, water_visc, time_step,
                              boundary_type='no_flow', flow_type='radial', width=None, length=None)
        df_tankcenter2['AWe'] = df_tankcenter2.get("AWe", w['Delta We'].to_list())
        df_tankcenter2['Cwe'] = df_tankcenter2.get("Cwe", w['Cumulative We'].to_list())
        we = df_tankcenter2['Cwe']
        EBM(p,np,wp,bo,cf,cw,sw0,boi, name,we,date)
    else:
        print("Por favor digita una opcion correcta")

    menu = int(input(
        "Menu: 1.tank_north(avg) \n 2. tank_center \n 3.tank_center(avg)"))

