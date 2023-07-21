import matplotlib.ticker as ticker
import streamlit as st
from streamlit_option_menu import option_menu
import pandas as pd
from pytank.aquifer import aquifer_carter_tracy,aquifer_fetkovich
from explore.test1 import Campbell,G_method
from pytank.utilities import days_in_month
from pytank.utilities import interp_from_dates, interp_dates_row
import matplotlib.pyplot as plt
formatter = ticker.EngFormatter()
st.title(" Material balance application ")
st.write(" This application ... ")

st.sidebar.title("DASHBOARD")
file=st.sidebar.file_uploader("Upload your file csv")

if file is None:
    st.write("Escoje un archivo")
else:
    df = pd.read_csv(file)
    df["Date"] = pd.to_datetime(df["Date"])
    st.write(df)
with st.sidebar:
    options = option_menu(
        menu_title="Main Menu",
        options=["Home","Exploratory Data Analysis" ,"Aquifer Calculation", "Graphic Material Balance"]
        ,
    )

if options == "Exploratory Data Analysis":
    file_prod = st.file_uploader("Upload you Production file csv")
    file_pvt = st.file_uploader("Upload your Pvt file csv")
    file_press = st.file_uploader("Upload your Pressure file csv")
    df_prod = pd.read_csv(file_prod)
    df_pvt = pd.read_csv(file_pvt)
    df_press = pd.read_csv(file_press)
    date_col = "START_DATETIME"
    df_prod[date_col] = pd.to_datetime(df_prod['START_DATETIME'])

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

    fig_1, (ax_11, ax_12) = plt.subplots(2, 1, sharex=True)

    (df_prod.pivot_table(oil_rate_col, date_col, well_name_col)
     .plot(colormap="Greens_r", lw=1, ax=ax_11, legend=False))
    (df_prod.pivot_table(water_rate_col, date_col, well_name_col)
     .plot(colormap="Blues_r", lw=1, ax=ax_12, legend=False))

    fig_1.suptitle("Production Rate per Well")
    ax_11.set_ylabel("Oil Rate (STB/D)")
    ax_12.set_ylabel("Water Rate (STB/D)")
    ax_12.set_xlabel("Date")
    st.pyplot(fig_1)

    fig_2, (ax_21, ax_22) = plt.subplots(2, 1, sharex=True)

    df_prod_tank = (df_prod[[date_col, tank_name_col, *cols_output]]
                    .groupby([date_col, tank_name_col])
                    .sum().reset_index())

    df_prod_tank.pivot_table(oil_rate_col, date_col, tank_name_col).plot(lw=1, ax=ax_21)
    df_prod_tank.pivot_table(water_rate_col, date_col, tank_name_col).plot(lw=1, ax=ax_22,
                                                                           legend=False)

    ax_21.legend(fontsize=6)
    fig_2.suptitle("Production Rate per Tank")
    ax_21.set_ylabel("Oil Rate (STB/D)")
    ax_22.set_ylabel("Water Rate (STB/D)")
    ax_22.set_xlabel("Date")
    st.pyplot(fig_2)

    # Rename column names for pressure data frame to use the same as the production one
    df_press.rename(columns={"WELLBORE": well_name_col, "DATE": date_col}, inplace=True)
    # Make sure the date column is o datetime object
    df_press[date_col] = pd.to_datetime(df_press[date_col])
    # Specify important columns
    press_col = "PRESSURE_DATUM"
    press_type_col = "TEST_TYPE"

    # For that, first calculate the liquid volume per well in each month
    liquid_vol_col = "liquid_vol"
    df_prod[liquid_vol_col] = df_prod[liquid_rate_col] * df_prod[cal_day_col]
    # Then, group by dates and sum the monthly volumes, then accumulate them
    df_field = df_prod.groupby(date_col)[liquid_vol_col].sum().cumsum().reset_index()
    # The resulting column is the liquid cumulative, so rename it
    df_field.rename(columns={liquid_vol_col: liquid_cum_col}, inplace=True)

    df_press[liquid_cum_col] = interp_from_dates(df_press[date_col],
                                                 df_field[date_col],
                                                 df_field[liquid_cum_col])

    fig_3, (ax_31, ax_32) = plt.subplots(2, 1)

    df_press.pivot_table(press_col, date_col, press_type_col).plot(style="o", ax=ax_31, ms=2)
    df_press.pivot_table(press_col, liquid_cum_col, press_type_col).plot(style="o", ax=ax_32,
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
    st.pyplot(fig_3)

    cols_group = [date_col, tank_name_col, liquid_vol_col]

    df_tank: pd.DataFrame = (
        df_prod[cols_group]
        .groupby(cols_group[:-1]).sum()  # to get monthly volumes per tank
        .groupby(tank_name_col).cumsum()  # to get cumulative vols per tank
        .reset_index())

    # The resulting column is the liquid cumulative, so rename it
    df_tank.rename(columns={liquid_vol_col: liquid_cum_col}, inplace=True)
    # Plot cumulatives to check
    ax_4 = (df_tank.pivot_table(liquid_cum_col, date_col, tank_name_col)
            .fillna(method="ffill")
            .plot())

    ax_4.set_title("Liquid Cumulatives per Tank")
    ax_4.set_ylabel("Liquid Cum (STB/D)")
    ax_4.set_xlabel("Date")
    ax_4.yaxis.set_major_formatter(formatter)
    #st.pyplot(ax_4)

    liquid_cum_col_per_tank = liquid_cum_col + "_tank"
    df_press[liquid_cum_col_per_tank] = (
        df_press.apply(lambda g: interp_dates_row(g, date_col, df_tank, date_col,
                                                  liquid_cum_col, tank_name_col,
                                                  tank_name_col), axis=1)
    )

    fig_5, (ax_51, ax_52) = plt.subplots(2, 1)

    df_press.pivot_table(press_col, date_col, tank_name_col).plot(ax=ax_51, style="o")
    df_press.pivot_table(press_col,
                         liquid_cum_col_per_tank, tank_name_col).plot(ax=ax_52,
                                                                      style="o",
                                                                      legend=False)

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
    st.pyplot(fig_5)

elif options == "Aquifer Calculation":
    select = st.selectbox("Aquifer Method",["Fetkovich","Carter-Tracy"])
    if select == "Fetkovich":
        aq_radius = int(st.number_input("Enter the value of the radius of the aquifer"))
        res_radius = int(st.number_input("Enter the value of the radius of the reservoir"))
        aq_thickness = int(st.number_input("Enter the value of the thickness of the aquifer"))
        phi = st.number_input("Enter the porosity value")
        theta = int(st.number_input("Enter the angle value"))
        k = int(st.number_input("Enter the permeability value"))
        ct = st.number_input("Enter the value of the total compressibility", format="%.6f")
        water_visc = st.number_input("Enter the value of the viscosity of the water")
        se1 = st.selectbox("Pressure",df.columns)
        pr=df[se1].values
        se2 = st.selectbox("Time step",df.columns)
        time_step = df[se2].values
        we=aquifer_fetkovich(aq_radius, res_radius, aq_thickness, phi, ct, pr, theta, k, water_visc, time_step,
            boundary_type='no_flow', flow_type='radial', width=None, length=None)
        df['Cwe'] = df.get("Cwe", we['Cumulative We'].to_list())
        st.write(we)
        st.write (df)
    elif select == "Carter-Tracy":
        res_radius = int(st.number_input("Enter the value of the radius of the reservoir"))
        aq_thickness = int(st.number_input("Enter the value of the thickness of the aquifer"))
        aq_por = st.number_input("Enter the porosity value")
        theta = int(st.number_input("Enter the angle value"))
        aq_perm = int(st.number_input("Enter the permeability value"))
        ct = st.number_input("Enter the value of the total compressibility", format="%.6f")
        water_visc = st.number_input("Enter the value of the viscosity of the water")
        se1 = st.selectbox("Pressure", df.columns)
        pr = df[se1].values
        se2 = st.selectbox("Time step", df.columns)
        time = df[se2].values
        we = aquifer_carter_tracy(aq_por, ct, res_radius, aq_thickness, theta, aq_perm,
                         water_visc, pr, time)
        df['Cwe'] = df.get("Cwe", we['Cumulative water influx, bbl'].to_list())
        st.write(we)

elif options== "Graphic Material Balance":
    st.write(df)
    df = df.drop(index=0)
    select = st.selectbox("Method", ["Campbell", "Graphical Method","Analytical method"])
    if select == "Campbell":
        se_p = st.selectbox("Pressure", df.columns)
        se_np = st.selectbox("Production oil cumulate", df.columns)
        se_wp = st.selectbox("Production water cumulate", df.columns)
        se_date = st.selectbox("Date", df.columns)
        se_bo = st.selectbox("Factor oil", df.columns)
        df[se_date]=df[se_date].dt.year
        p = df[se_p].values
        np = df[se_np].values
        wp = df[se_wp].values
        bo = df[se_bo].values
        date = df[se_date].values
        cf = st.number_input("Enter the value of the formatio compressibility", format="%.6f")
        cw = st.number_input("Enter the value of the water compressibility", format="%.6f")
        sw0 = st.number_input("Enter the value of the Saturation water")
        pi = st.number_input("Enter the value Inicial Pressure")
        boi = st.number_input("Enter the value of the total compressibility", format="%.4f")

        gra= Campbell(p, np,wp, bo, cf, cw, sw0, boi, date,pi)
        st.pyplot(gra)

    elif select== "Graphical Method":
        se_p = st.selectbox("Pressure", df.columns)
        se_np = st.selectbox("Production oil cumulate", df.columns)
        se_wp = st.selectbox("Production water cumulate", df.columns)
        se_we = st.selectbox("Cumulative We", df.columns)
        se_bo = st.selectbox("Factor oil", df.columns)
        pr = df[se_p].values
        np = df[se_np].values
        wp = df[se_wp].values
        bo = df[se_bo].values
        we = df[se_we].values
        cf = st.number_input("Enter the value of the formatio compressibility", format="%.6f")
        cw = st.number_input("Enter the value of the water compressibility", format="%.6f")
        sw0 = st.number_input("Enter the value of the Saturation water")
        pi = st.number_input("Enter the value Inicial Pressure")
        boi = st.number_input("Enter the value Inicial Bo", format="%.4f")

        gra= G_method(pr, np,wp, bo, cf, cw, sw0, boi,pi,we)
        st.pyplot(gra)