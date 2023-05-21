#%%
from pytank.pvt_table import pvt_table
import matplotlib.pyplot as plt

#%%
p_sep = 80  # Separator pressure [psia]
t_sep = 80  # Separator temperature [F]
api = 38  # Stock tank oil gravity [API]
rsp = 520  # Separator producing GOR [scf/STB]
sg_sep = 0.8  # Separator gas specific gravity
tres = 210  # Reservoir temperature [F]
den_sto = 40  # Density of stock-tank oil [lb/cu ft]
p_res = 3000  # Initial reservoir pressure [psia]
salinity = 35000  # Water Salinity [ppm]
Jump = 100
units=1
#%%
tabla = pvt_table(p_sep, t_sep, api, rsp, sg_sep, tres, den_sto, p_res, salinity, Jump,
                  units)
