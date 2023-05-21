#%%
from pytank.material_balance import campbell_function, ho_terms_equation, \
    havlena_and_odeh
import pandas as pd


#%% Load the the mbal dataframe
df_mbal = pd.read_csv("mbal_Dataframe.csv")

#%% Creation of a mbal dataframe for each tank
df_tankcenter = df_mbal[df_mbal["Tank"] == "tank_center"]
df_tanksouth = df_mbal[df_mbal["Tank"] == "tank_south"]
df_tanknorth = df_mbal[df_mbal["Tank"] == "tank_north"]

#%% Additional numerical data
water_sat = 0.15
rock_comp = 3.5e-6
water_comp = 3.62e-6
water_fvf = 1.0
gas_water_rs = 0
oil_fvf_init = 0.86
gas_fvf_init = 0.04
tot_fvf_init = 0.60
gas_oil_rs_init = 89
pressure_init = 3000

#%% Calculation of the Havlena and Odeh terms for each tank, as well as its respective
# Campbell plot

""" ----------------------------TANK NORTH-------------------------------------------"""

df_mbal_tankNorth = ho_terms_equation(df_tanknorth, "oil_prod_cum", "water_prod_cum",
                             "gas_prod_cum", "PRESSURE_DATUM", "oil_fvf", "gas_fvf",
                             "gas_oil_rs", water_fvf, gas_water_rs,
                        water_sat, water_comp, rock_comp, oil_fvf_init, gas_fvf_init,
                                      tot_fvf_init, gas_oil_rs_init, pressure_init)

#%%
campbell_plot_tankNorth = campbell_function(df_tanknorth, "oil_prod_cum",
                                             "water_prod_cum", "gas_prod_cum",
                                             "PRESSURE_DATUM", "F", "Eo", "Efw",
                                             "oil_fvf", "gas_fvf", "gas_oil_rs",
                                             water_fvf, gas_water_rs, water_sat,
                                             water_comp, rock_comp, oil_fvf_init,
                                            gas_fvf_init, tot_fvf_init, gas_oil_rs_init,
                                            pressure_init)


#%%
""" ----------------------------TANK CENTER------------------------------------------"""

df_mbal_tankCenter = ho_terms_equation(df_tankcenter, "oil_prod_cum", "water_prod_cum",
                             "gas_prod_cum", "PRESSURE_DATUM", "oil_fvf", "gas_fvf",
                             "gas_oil_rs", water_fvf, gas_water_rs,
                        water_sat, water_comp, rock_comp, oil_fvf_init, gas_fvf_init,
                                       tot_fvf_init, gas_oil_rs_init, pressure_init)

#%%
campbell_plot_tankCenter = campbell_function(df_tankcenter, "oil_prod_cum",
                                             "water_prod_cum", "gas_prod_cum",
                                             "PRESSURE_DATUM", "F", "Eo", "Efw",
                                             "oil_fvf", "gas_fvf", "gas_oil_rs",
                                             water_fvf, gas_water_rs, water_sat,
                                             water_comp, rock_comp, oil_fvf_init,
                                             gas_fvf_init, tot_fvf_init, gas_oil_rs_init,
                                             pressure_init)


#%%
""" ----------------------------TANK SOUTH-------------------------------------------"""

df_mbal_tankSouth = ho_terms_equation(df_tanksouth, "oil_prod_cum", "water_prod_cum",
                             "gas_prod_cum", "PRESSURE_DATUM", "oil_fvf", "gas_fvf",
                             "gas_oil_rs", water_fvf, gas_water_rs,
                        water_sat, water_comp, rock_comp, oil_fvf_init, gas_fvf_init,
                                      tot_fvf_init, gas_oil_rs_init, pressure_init)

#%%
campbell_plot_tankSouth = campbell_function(df_tanksouth, "oil_prod_cum",
                                             "water_prod_cum", "gas_prod_cum",
                                             "PRESSURE_DATUM", "F", "Eo", "Efw",
                                             "oil_fvf", "gas_fvf", "gas_oil_rs",
                                             water_fvf, gas_water_rs, water_sat,
                                             water_comp, rock_comp, oil_fvf_init,
                                            gas_fvf_init, tot_fvf_init, gas_oil_rs_init,
                                            pressure_init)





