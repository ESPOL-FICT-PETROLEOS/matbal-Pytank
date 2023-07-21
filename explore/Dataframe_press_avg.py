#%%
import pandas as pd
from pytank.utilities import interp_dates_row
from pytank.pvt_interp import interp_pvt_matbal
import matplotlib.pyplot as plt
from pytank.material_balance import underground_withdrawal

# Load data into dataframes
df_pressure = pd.read_csv("../tests/data_for_tests/full_example_1/pressures.csv")
df_prod = pd.read_csv("../tests/data_for_tests/full_example_1/production.csv")
df_avg=pd.read_csv("press_avg.csv")
df_avg['PRESSURE_DATUM']=df_avg['PRESSURE_DATUM'].interpolate(method="linear")
#%% Casting of date data type of pressure and prod dataframes
date_col = "DATE"
datep_col = "START_DATETIME"
df_avg[date_col] = pd.to_datetime(df_avg['START_DATETIME'])
df_prod[datep_col] = pd.to_datetime(df_prod['START_DATETIME'])

# %% Define data frame columns of production dataframe
# Input
oil_cum_col = "OIL_CUM"
water_cum_col = "WATER_CUM"
gas_cum_col = "GAS_CUM"
well_name_col = "ITEM_NAME"
tank_name_col = "Tank"
# output
oil_vol_col = "oil_vol"
water_vol_col = "water_col"
gas_vol_col = "gas_col"

#%% Calculate monthly volume
cols_input = [oil_cum_col, water_cum_col, gas_cum_col]
cols_output = [oil_vol_col, water_vol_col, gas_vol_col]

df_prod[cols_output] = (df_prod.groupby(well_name_col)[cols_input]
                        .diff().fillna(df_prod[cols_input]))

#%% Creation of a new dataframe grouping the cumulative fluids produced by each tank
cols_group = [datep_col, tank_name_col, oil_vol_col, water_vol_col, gas_vol_col]

df_tanks = (
    df_prod[cols_group]
    .groupby(cols_group[0:2])
    .sum()
    .groupby(tank_name_col)
    .cumsum()
    .reset_index()
)

df_tanks.rename(columns={oil_vol_col: oil_cum_col,
                         water_vol_col: water_cum_col,
                         gas_vol_col: gas_cum_col}, inplace=True)

#%% Plot Oil Cumulative
df_tanks.pivot_table(oil_cum_col,datep_col, tank_name_col).fillna(method="ffill").plot()
plt.show()
#df_tanks.pivot(datep_col, tank_name_col, oil_cum_col).fillna(method="ffill").plot()
#plt.show()

# %% Housekeeping of pressure data frame
# Rename column names for pressure data frame to use the same as the production one
df_avg.rename(columns={"WELLBORE": well_name_col, "DATE": date_col}, inplace=True)

#%% Using the interp_dates_row function to add a new column to the pressure Dataframe
# named oil_cum_col per tank
oil_cum_col_per_tank = oil_cum_col + "_Tank"
df_avg[oil_cum_col_per_tank] = df_avg.apply(
    lambda g: interp_dates_row(
        g, date_col, df_tanks, datep_col, oil_cum_col, tank_name_col, tank_name_col
    ),
    axis=1,
)

#%% Using the interp_dates_row function to add a new column to the pressure Dataframe
# named water_cum_col per tank
water_cum_col_per_tank = water_cum_col + "_Tank"
df_avg[water_cum_col_per_tank] = df_avg.apply(
    lambda g: interp_dates_row(
        g, date_col, df_tanks, datep_col, water_cum_col, tank_name_col, tank_name_col
    ),
    axis=1,
)

#%% Using the interp_dates_row function to add a new column to the pressure Dataframe
# named gas_cum_col per tank
gas_cum_col_per_tank = gas_cum_col + "_Tank"
df_avg[gas_cum_col_per_tank] = df_avg.apply(
    lambda g: interp_dates_row(
        g, date_col, df_tanks, datep_col, gas_cum_col, tank_name_col, tank_name_col
    ),
    axis=1,
)

#%% Deletion of unnecessary columns, as well as, sorting by time, and finally, the
# pressure dataframe is renamed as mbal

Sand = "SAND"
Test_type = "TEST_TYPE"

df_mbal2 = df_avg
df_mbal2 = df_mbal2.sort_values(date_col)

#%% Creation of a mbal dataframe for each tank
df_center = df_mbal2[df_mbal2[tank_name_col] == "tank_center"]
df_south = df_mbal2[df_mbal2[tank_name_col] == "tank_south"]
df_north = df_mbal2[df_mbal2[tank_name_col] == "tank_north"]

#%% Load of the pvt dataframe, which is the one used to interpolate and then, add
# the pvt columns in the mbal dataframe
df_pvt = pd.read_csv("../tests/data_for_tests/full_example_1/pvt.csv")

# filling the missing values of the Bg column using the ffill method
df_pvt = df_pvt.fillna(method="ffill")

# Define the names of the new columns to be added in the mbal dataframe, which
# corresponds to the pvt properties
oil_fvf_col = "Bo"
gas_fvf_col = "Bg"
gas_oil_rs_col = "GOR"
ppvt_col = "Pressure"
pres_col = "PRESSURE_DATUM"

# Creation of empty lists for each pvt column
oil_fvf = []
gas_fvf = []
gas_o_rs = []

#%% Loop into the pressure column of the mbal dataframe in order to use the
# interp_pvt_matbal function to interpolate and calculate the oil_fvf corresponding to
# each of its pressure values

for press in df_center[pres_col]:
    oil_fvf.append(interp_pvt_matbal(df_pvt, ppvt_col, oil_fvf_col, press))

# Transforming this numpy array to a pandas Series
oil_fvf = pd.Series(oil_fvf, name="oil_fvf")


#%% # Loop in the pressure column of the mbal dataframe in order to use the
# interp_pvt_matbal function to interpolate and calculate the gas_fvf corresponding to
# each of its pressure values

for press in df_center[pres_col]:
    gas_fvf.append(interp_pvt_matbal(df_pvt, ppvt_col, gas_fvf_col, press))

# Transforming this numpy array to a pandas Series
gas_fvf = pd.Series(gas_fvf, name="gas_fvf")

#%% # Loop in the pressure column of the mbal dataframe in order to use the
# interp_pvt_matbal function to interpolate and calculate the gas_oil_rs_col corresponding
# to each of its pressure values

for press in df_center[pres_col]:
    gas_o_rs.append(interp_pvt_matbal(df_pvt, ppvt_col, gas_oil_rs_col, press))

# Transforming this numpy array to a pandas Series
gas_oil_rs = pd.Series(gas_o_rs, name="gas_oil_rs_col")

#%% Concatenation of the mbal dataframe with its respective pvt columns (Pandas Series)
df_center = pd.concat([df_center, oil_fvf, gas_fvf, gas_oil_rs], axis=1).sort_values(
    date_col
)

#%% Rename of the columns of the mbal dataframe to match the variables names stated
# in the file attached in unterbase
df_center.rename(
    columns={
        "OIL_CUM_Tank": "oil_prod_cum",
        "WATER_CUM_Tank": "water_prod_cum",
        "GAS_CUM_Tank": "gas_prod_cum","START_DATETIME": "Date","PRESSURE_DATUM": "Pressure"
    },
    inplace=True,
)
df_center =df_center.drop(columns=['DATE'])
nueva_fila = pd.DataFrame({'Date': '1987-09-01', 'Tank': 'tank_center', 'Pressure': 3700.00, 'oil_fvf': 1.1},
                          index=[0])
df_center = pd.concat([nueva_fila, df_center]).reset_index(drop=True)
df_center.iloc[0] = df_center.iloc[0].fillna(0)
df_center['Date']=pd.to_datetime(df_center['Date'])
fecha= pd.to_datetime("1987-09-01")
df_center['time_setp'] = (df_center['Date']-fecha).dt.days
#%% Convert the mbal dataframe to a csv file, which then will be imported to calculate
# the terms F, Eo, Eg, and Efw
mbal = df_center.to_csv("mbal_Dataframe3.csv", index=False)


#%% Using the interp_dates_row function to add the pressure column to the df_tanks
# Dataframe
Pressure = "Pressure"
df_tanks[Pressure] = df_tanks.apply(
    lambda g: interp_dates_row(
        g, datep_col, df_mbal2, date_col, pres_col, tank_name_col, tank_name_col
    ),
    axis=1,
)

#%% Loop into the pressure column of the df_tanks dataframe in order to use the
# interp_pvt_matbal function to interpolate and calculate the oil_fvf corresponding to
# each of its pressure values

Oil_fvf = []
for press in df_tanks[Pressure]:
    Oil_fvf.append(interp_pvt_matbal(df_pvt, ppvt_col, oil_fvf_col, press))

# Transforming this numpy array to a pandas Series
Oil_fvf = pd.Series(Oil_fvf, name="oil_fvf")


#%% # Loop in the pressure column of the df_tanks dataframe in order to use the
# interp_pvt_matbal function to interpolate and calculate the gas_fvf corresponding to
# each of its pressure values

Gas_fvf = []
for press in df_tanks[Pressure]:
    Gas_fvf.append(interp_pvt_matbal(df_pvt, ppvt_col, gas_fvf_col, press))

# Transforming this numpy array to a pandas Series
Gas_fvf = pd.Series(Gas_fvf, name="gas_fvf")


#%% # Loop in the pressure column of the df_tanks dataframe in order to use the
# interp_pvt_matbal function to interpolate and calculate the gas_oil_rs_col corresponding
# to each of its pressure values

Gas_o_rs = []
for press in df_tanks[Pressure]:
    Gas_o_rs.append(interp_pvt_matbal(df_pvt, ppvt_col, gas_oil_rs_col, press))

# Transforming this numpy array to a pandas Series
Gas_o_rs = pd.Series(Gas_o_rs, name="gas_oil_rs_col")


#%% Concatenation of the mbal dataframe with its respective pvt columns (Pandas Series)
df_tanks = pd.concat([df_tanks, Oil_fvf, Gas_fvf, Gas_o_rs], axis=1)


#%% Rename of the columns of the df_tanks dataframe to match the variables names stated
# in the file attached in unterbase
df_tanks.rename(
    columns={
        "OIL_CUM": "oil_prod_cum",
        "WATER_CUM": "water_prod_cum",
        "GAS_CUM": "gas_prod_cum",
    "START_DATETIME": "Date"},
    inplace=True,
)

#%% Convert the df_tanks dataframe to a csv file, which then will be imported to
# calculate the terms F, Eo, Eg, and Efw
mbal_tanks = df_tanks.to_csv("mbal_tanks2.csv", index=False)

#%% New Dataframe with pressure avg

