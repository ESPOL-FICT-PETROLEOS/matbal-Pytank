import pandas as pd
import matplotlib.pyplot as plt
from pytank.utilities import days_in_month
from pytank.utilities import interp_from_dates, interp_dates_row
import matplotlib.ticker as ticker
from scipy.interpolate import interp1d
from pytank.material_balance import underground_withdrawal, pressure_vol_avg

formatter = ticker.EngFormatter()

# %% specify files to load
production_file = "tests/data_for_tests/full_example_1/production.csv"
pressure_file = "tests/data_for_tests/full_example_1/pressures.csv"
pvt_file = "tests/data_for_tests/full_example_1/pvt.csv"

# %% Load data into dataframes
df_prod = pd.read_csv(production_file)
df_press = pd.read_csv(pressure_file)
df_pvt = pd.read_csv(pvt_file)

# %% Cast date column to date
date_col = "START_DATETIME"
df_prod[date_col] = df_prod[date_col].astype("datetime64")

# %% Define data frame columns
# Input
oil_cum_col = "OIL_CUM"
water_cum_col = "WATER_CUM"
gas_cum_col = "GAS_CUM"
well_name_col = "ITEM_NAME"
tank_name_col = "Tank"
# Output
cal_day_col = "cal_day"
oil_rate_col = "oil_rate"
water_rate_col = "water_rate"
gas_rate_col = "gas_rate"
liquid_rate_col = "liquid_rate"
liquid_cum_col = "liquid_cum"

# %% Calculate Rates
# Calculate the calendar days
df_prod[cal_day_col] = df_prod[date_col].map(lambda date: days_in_month(date))

# Define the input and output columns
cols_input = [oil_cum_col, water_cum_col, gas_cum_col]
cols_output = [oil_rate_col, water_rate_col, gas_rate_col]

# Calculate the rates using the differences between cumulatives
df_input = df_prod[[well_name_col, *cols_input]]
df_prod[cols_output] = (df_input.groupby(well_name_col).diff().fillna(df_input)
                        .div(df_prod[cal_day_col], axis=0))

# Calculate liquid production rate and cumulative
df_prod[liquid_rate_col] = df_prod[oil_rate_col] + df_prod[water_rate_col]
df_prod[liquid_cum_col] = df_prod[oil_cum_col] + df_prod[water_cum_col]

# %% Plot oil production from all the wells
fig_1, (ax_11, ax_12) = plt.subplots(2, 1, sharex=True)

(df_prod.pivot(date_col, well_name_col, oil_rate_col)
 .plot(colormap="Greens_r", lw=1, ax=ax_11, legend=False))
(df_prod.pivot(date_col, well_name_col, water_rate_col)
 .plot(colormap="Blues_r", lw=1, ax=ax_12, legend=False))

fig_1.suptitle("Production Rate per Well")
ax_11.set_ylabel("Oil Rate (STB/D)")
ax_12.set_ylabel("Water Rate (STB/D)")
ax_12.set_xlabel("Date")
plt.show()

# %% Plot the Oil and Liquid Production per tank
fig_2, (ax_21, ax_22) = plt.subplots(2, 1, sharex=True)

df_prod_tank = (df_prod[[date_col, tank_name_col, *cols_output]]
                .groupby([date_col, tank_name_col])
                .sum().reset_index())

df_prod_tank.pivot(date_col, tank_name_col, oil_rate_col).plot(lw=1, ax=ax_21)
df_prod_tank.pivot(date_col, tank_name_col, water_rate_col).plot(lw=1, ax=ax_22,
                                                                 legend=False)

ax_21.legend(fontsize=6)
fig_2.suptitle("Production Rate per Tank")
ax_21.set_ylabel("Oil Rate (STB/D)")
ax_22.set_ylabel("Water Rate (STB/D)")
ax_22.set_xlabel("Date")
plt.show()

# %% Housekeeping of pressure data frame
# Rename column names for pressure data frame to use the same as the production one
df_press.rename(columns={"WELLBORE": well_name_col, "DATE": date_col}, inplace=True)
# Make sure the date column is o datetime object
df_press[date_col] = df_press[date_col].astype("datetime64")
# Specify important columns
press_col = "PRESSURE_DATUM"
press_type_col = "TEST_TYPE"

# %% Calculate liquid cumulative in the whole field
# For that, first calculate the liquid volume per well in each month
liquid_vol_col = "liquid_vol"
df_prod[liquid_vol_col] = df_prod[liquid_rate_col] * df_prod[cal_day_col]
# Then, group by dates and sum the monthly volumes, then accumulate them
df_field = df_prod.groupby(date_col)[liquid_vol_col].sum().cumsum().reset_index()
# The resulting column is the liquid cumulative, so rename it
df_field.rename(columns={liquid_vol_col: liquid_cum_col}, inplace=True)

# %% Prepare pressure data for plotting
df_press[liquid_cum_col] = interp_from_dates(df_press[date_col],
                                             df_field[date_col],
                                             df_field[liquid_cum_col])

# %% Plot all pressure data colouring by type, vs date and vs Liq. Cum.
fig_3, (ax_31, ax_32) = plt.subplots(2, 1)

df_press.pivot(date_col, press_type_col, press_col).plot(style="o", ax=ax_31, ms=2)
df_press.pivot(liquid_cum_col, press_type_col, press_col).plot(style="o", ax=ax_32,
                                                               ms=2, legend=False)
# Set first axes
ax_31.set_title("Pressure data vs. Date")
ax_31.set_xlabel("Date")
ax_31.set_ylabel("Pressure (psia)")
ax_31.tick_params(axis="x", labelsize=8)
ax_31.legend(fontsize=8)
ax_31.yaxis.set_major_formatter(formatter)
# Set second axes
ax_32.set_title("Pressure data vs. Liquid Cumulative")
ax_32.set_xlabel("Liquid Cumulative (STB)")
ax_32.set_ylabel("Pressure (psia)")
ax_32.xaxis.set_major_formatter(formatter)
ax_32.yaxis.set_major_formatter(formatter)

plt.tight_layout()
plt.show()

# %% Calculate liquid cumulative per tank in the pressure data frame
cols_group = [date_col, tank_name_col, liquid_vol_col]

df_tank: pd.DataFrame = (
    df_prod[cols_group]
        .groupby(cols_group[:-1]).sum()  # to get monthly volumes per tank
        .groupby(tank_name_col).cumsum()  # to get cumulative vols per tank
        .reset_index())

# The resulting column is the liquid cumulative, so rename it
df_tank.rename(columns={liquid_vol_col: liquid_cum_col}, inplace=True)
# Plot cumulatives to check
ax_4 = (df_tank.pivot(date_col, tank_name_col, liquid_cum_col)
        .fillna(method="ffill")
        .plot())

ax_4.set_title("Liquid Cumulatives per Tank")
ax_4.set_ylabel("Liquid Cum (STB/D)")
ax_4.set_xlabel("Date")
ax_4.yaxis.set_major_formatter(formatter)
plt.show()

# %% Get liquid cumulative interpolation per tank into the pressure data frame
liquid_cum_col_per_tank = liquid_cum_col + "_tank"
df_press[liquid_cum_col_per_tank] = (
    df_press.apply(lambda g: interp_dates_row(g, date_col, df_tank, date_col,
                                              liquid_cum_col, tank_name_col,
                                              tank_name_col), axis=1)
)

#%% Plot pressure per tank vs. date and vs. liquid cumulative
fig_5, (ax_51, ax_52) = plt.subplots(2, 1)

df_press.pivot_table(press_col, date_col, tank_name_col).plot(ax=ax_51, style="o")
df_press.pivot_table(press_col,
                     liquid_cum_col_per_tank,
                     tank_name_col).plot(ax=ax_52,
                                         style="o",
                                         legend=False)

# Set first axes
ax_51.set_title("Pressure data vs. Date")
ax_51.set_xlabel("Date")
ax_51.set_ylabel("Pressure (psia)")
ax_51.tick_params(axis="x", labelsize=8)
ax_51.legend(fontsize=8)
ax_51.yaxis.set_major_formatter(formatter)
# Set second axes
ax_52.set_title("Pressure data vs. Liquid Cumulative")
ax_52.set_xlabel("Liquid Cumulative (STB)")
ax_52.set_ylabel("Pressure (psia)")
ax_52.xaxis.set_major_formatter(formatter)
ax_52.yaxis.set_major_formatter(formatter)

plt.tight_layout()
plt.show()

#%% ------------- Calculate pressure volumetric average --------------------------------


#%% 1.- Interpolate PVT information

# Define the fluid properties column names for the df pressure
oil_fvf_col = "Bo"
gas_oil_rs_col = "Rs"
gas_fvf_col = "Bg"

oil_fvf_interp = interp1d(df_pvt["Pressure"],
                          df_pvt["Bo"],
                          fill_value="extrapolate")

gas_oil_rs_interp = interp1d(df_pvt["Pressure"],
                             df_pvt["GOR"],
                             fill_value="extrapolate")

gas_fvf_interp = interp1d(df_pvt["Pressure"],
                          df_pvt["Bg"],
                          fill_value="extrapolate")

# 2.- Apply the functions to the pressure df
df_press[oil_fvf_col] = oil_fvf_interp(df_press[press_col])
df_press[gas_oil_rs_col] = gas_oil_rs_interp(df_press[press_col])
df_press[gas_fvf_col] = gas_fvf_interp(df_press[press_col])

#%% Interpolate oil, gas and water cumulatives into the pressure data frame

for col in cols_input:
    df_press[col] = df_press.apply(
        lambda x: interp_dates_row(x, date_col, df_prod, date_col,
                                   col, well_name_col, well_name_col,
                                   left=0.0), axis=1
    )
    # For wells not available in the production data frame, fill nans with 0
    df_press[col].fillna(0, inplace=True)

#%% Calculate underground withdrawal for each well
uw_col = "UW"

df_press[uw_col] = underground_withdrawal(df_press, oil_cum_col, water_cum_col,
                                          gas_cum_col, oil_fvf_col, 1,
                                          gas_fvf_col, gas_oil_rs_col, 0)

#%% Calculate the pressure volumetric average per tank
avg_freq = "12MS"
df_press_avg = df_press.groupby(tank_name_col).apply(
    lambda g: pressure_vol_avg(g, well_name_col, date_col, press_col, uw_col,
                               avg_freq, "end")
).reset_index(0)
#%% Plot the pressure average and normal pressures of the central tank
df_press_avg_tank = df_press_avg.loc[df_press_avg[tank_name_col] == "tank_center",
                                     [date_col, press_col]]

df_press_tank = df_press.loc[df_press[tank_name_col] == "tank_center",
                             [date_col, press_col]]

fig_6, ax_6 = plt.subplots(1, 1)

df_press_avg_tank.plot(x=date_col, y=press_col, ax=ax_6, style="bo", label="avg")
df_press_tank.plot(x=date_col, y=press_col, ax=ax_6, style="ro", label="data")
plt.show()
