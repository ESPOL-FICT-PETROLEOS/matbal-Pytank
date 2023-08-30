import matplotlib.ticker as ticker
import streamlit as st
from streamlit_option_menu import option_menu
import pandas as pd
from pytank.aquifer import aquifer_carter_tracy,aquifer_fetkovich
from explore.funciones_app import Campbell,G_method,calcuted_pressure,cw,bw
from pytank.utilities import days_in_month
from pytank.utilities import interp_from_dates, interp_dates_row
import matplotlib.pyplot as plt
import plotly.express as px
from scipy import stats
from PIL import Image
formatter = ticker.EngFormatter()
icon = Image.open("resources/icon.PNG")
st.set_page_config(page_title="Material balance",page_icon=icon)
st.title(" Material balance application :memo:")
logo=Image.open("resources/logo_nuevo.png")
st.sidebar.image(logo)
st.sidebar.title("Dashboard :large_orange_diamond:")

try:
    file = st.sidebar.file_uploader("Upload your Dataframe file csv")
    file_pvt = st.sidebar.file_uploader("Upload your Pvt file csv")

    if file is not None and file_pvt is not None:
        df = pd.read_csv(file)
        df_pvt = pd.read_csv(file_pvt)
        df["Date"] = pd.to_datetime(df["Date"])
        st.write(df)

except Exception as e:
    st.error(f"An error occurred: {e}")

with st.sidebar:
    options = option_menu(
        menu_title="Main Menu",
        options=["Home","Exploratory Data Analysis" ,"Graphic Campbell", "Graphical Method","Analytical method"]
        , icons=["house-gear","bar-chart-fill","bar-chart-steps","graph-up","graph-down-arrow"],
    )
if options == "Home":
    st.header("Material Balance Analysis Theory")
    image5 = Image.open("resources/oil reservoir.png")
    st.image(image5)
    st.write("Material balance analysis is an interpretation method used to determine original fluids-in-place (OFIP) based on production and static pressure data."
             "The general material balance equation relates the original oil, gas, and water in the reservoir to production volumes and current pressure conditions / fluid properties. The material balance equations considered assume tank type behaviour at any given datum depth - the reservoir is considered to have the same pressure and fluid properties at any location in the reservoir. This assumption is quite reasonable provided that quality production and static pressure measurements are obtained."
             "Consider the case of the depletion of the reservoir pictured below. At a given time after the production of fluids from the reservoir has commenced, the pressure will have dropped from its initial reservoir pressure pi, to some average reservoir pressure p. Using the law of mass balance, during the pressure drop (Dp), the expansion of the fluids leftover in the reservoir must be equal to the volume of fluids produced from the reservoir.")
    st.write("The simplest way to visualize material balance is that if the measured surface volume of oil, gas and water were returned to a reservoir at the reduced pressure, it must fit exactly into the volume of the total fluid expansion plus the fluid influx."
             "The general form of the equation can be described as net withdrawal (withdrawal - injection) = expansion of the hydrocarbon fluids in the system + cumulative water influx. This is shown in the equation below:")
    latex_equation = r'''
    N_p \left[ B_o + (R_p - R_s) B_g \right] + W_p B_w = N \left[ (B_o - B_{oi}) + (R_{si} - R_s) B_g \right] + mN B_{oi} \left( \frac{B_g}{B_{gi}} - 1 \right) 
    + N B_{oi} \left( 1 + m \right) \left( \frac{c_w S_{wi} + c_f}{1 - S_{wi}} \right) \Delta p + W_e 
    '''
    st.latex(latex_equation)

    st.subheader("Diagnostics - Dake and Campbell")
    st.write("Dake and Campbell plots are used as diagnostic tools to identify the reservoir type based on the signature of production and pressure behaviour. The plots are established based on the assumption of a volumetric reservoir, and deviation from this behaviour is used to indicate the reservoir type."
             "In the Dake plot, the simplest oil case of solution gas / depletion drive (no gas cap, no water drive) is used to determine the axes of the plot. The material balance equation is rearranged as shown below.")
    latex_equation2 = r'''
    \frac{F}{{E_o + E_{f,w}}} = N
    '''
    st.latex(latex_equation2)
    image1=Image.open("resources/campbell.png")
    st.image(image1)
    st.subheader("Havlena-Odeh")
    st.write("As seen in the general material balance equation there are many unknowns, and as a result finding an exact or unique solution can be difficult. "
             "However, using other techniques to help determine some variables (for example, m or original gas-in-place from volumetrics or seismic), the equation can be simplified to yield a more useful answer. Various plots are available to conduct an oil material balance rather than calculating an answer from individual measurements of reservoir pressure.")
    latex_equation3 = r'''
    \frac{F}{{E_o + E_{f,w}}} = N + \frac{{We \cdot Bw}}{{E_o + E_{f,w}}}
    '''

    st.latex(latex_equation3)
    image2=Image.open("resources/odeh.jpg")
    st.image(image2, width=400)
    st.subheader("Pressure History Match")
    st.write("""
The pressure match method employs an iterative procedure that uses the values of original oil-in-place, original gas-in-place, and W to calculate the reservoir pressure versus time. The pressure match is then plotted against the real measured static reservoir pressures and compared. This is by far the most robust and easily understood material balance technique for the following reasons:

1. Pressure and time are easily understood variables, and so sensitivity analysis can be conducted relatively easily.

2. Use of time enables the analyst to see directly the impact of:

    Changing withdrawal rates, especially shut-ins on reservoir pressure decline.

    Water drive and connected reservoirs on reservoir depletion, especially since these are both cumulative withdrawal and time-based processes.

3. A relatively simple, iterative process is used to achieve a unique solution wherein:

    Start with the simplest solution (oil and/or gas depletion only) and then proceed to more complex models only if demonstrated to be required


An example of a pressure history match is shown below:
""")
    image3 = Image.open("resources/pvst.PNG")
    st.image(image3, width=350)
elif options == "Exploratory Data Analysis":
    try:
        file_prod = st.file_uploader("Upload you Production file csv")
        file_press = st.file_uploader("Upload your Pressure file csv")
        df_prod = pd.read_csv(file_prod)
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
    except Exception as e:
        st.write("")


elif options== "Graphic Campbell":
    try:
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
        cf = st.sidebar.number_input("formatio compressibility", format="%.7f")
        sw0 = st.sidebar.number_input("Saturation water")
        pi = st.sidebar.number_input("Inicial Pressure")
        boi = st.sidebar.number_input("Inicial Factor oil", format="%.5f")
        t = st.sidebar.number_input("Temperature")
        salinity = st.sidebar.number_input("Salinity")
        gra= Campbell(p, np,wp, bo, cf, sw0, boi, date,pi,t,salinity)
        fig = px.scatter(gra, x='Date', y='F/Eo+Efw', title='Gráfico Campbell')
        fig.update_xaxes(title_text='Date (Year)')
        fig.update_yaxes(title_text='F/Eo+Efw')
        x= gra['Date']
        y= gra['F/Eo+Efw']
        slope, intercept, r, p, se = stats.linregress(x, y)
        yt = intercept + (x * slope)
        fig.add_scatter(x=x, y=yt, mode='lines', line=dict(color='green'),
                        name='Regression Lineal')
        # Mostrar el gráfico en Streamlit
        st.plotly_chart(fig)
    except Exception as e:
        st.write("")

elif options== "Graphical Method":
    try:
        select = st.radio("Aquifer Method", ["Fetkovich", "Carter-Tracy"])
        t = st.sidebar.number_input("Temperature")
        salinity = st.sidebar.number_input("Salinity")
        cf = st.sidebar.number_input("Formation Compressibility", format="%.6f")
        se_p = st.selectbox("Pressure", df.columns)
        pr = df[se_p].values
        Cw = cw(pr, t, salinity)['Cw']
        if select == "Fetkovich":
            aq_radius = st.sidebar.number_input("Aquifer Radius ")
            res_radius = st.sidebar.number_input("Reservoir Radius")
            aq_thickness = st.sidebar.number_input("Aquifer Thickness")
            phi = st.sidebar.number_input("Porosity")
            theta = st.sidebar.number_input("Angle")
            k = st.sidebar.number_input("Permeability")
            water_visc = st.sidebar.number_input("Water viscosity ")
            se2 = st.selectbox("Time step", df.columns)
            time_step = df[se2].values
            we = aquifer_fetkovich(aq_radius, res_radius, aq_thickness, phi, cf, pr, theta, k, water_visc, time_step,
                                    boundary_type='no_flow', flow_type='radial', width=None, length=None)
            df['Cwe'] = df.get("Cwe", we['Cumulative We'].to_list())
            st.write(we)
        elif select == "Carter-Tracy":
            res_radius = st.sidebar.number_input("Reservoir Radius")
            aq_thickness = st.sidebar.number_input("Aquifer Thickness")
            aq_por = st.sidebar.number_input("Porosity")
            theta = st.sidebar.number_input("Angle")
            aq_perm = st.sidebar.number_input("Permeability")
            water_visc = st.sidebar.number_input("Water viscosity")
            se2 = st.selectbox("Time step", df.columns)
            time = df[se2].values
            we = aquifer_carter_tracy(aq_por, cf, res_radius, aq_thickness, theta, aq_perm,
                                      water_visc, pr, time)
            df['Cwe'] = df.get("Cwe", we['Cumulative water influx, bbl'].to_list())
            st.write(we)
        se_np = st.selectbox("Production oil cumulate", df.columns)
        se_wp = st.selectbox("Production water cumulate", df.columns)
        se_we = st.selectbox("Cumulative We", df.columns)
        se_bo = st.selectbox("Factor oil", df.columns)
        pr = df[se_p].values
        np = df[se_np].values
        wp = df[se_wp].values
        bo = df[se_bo].values
        we = df[se_we].values
        sw0 = st.sidebar.number_input("Saturation water")
        pi = st.sidebar.number_input("Inicial Pressure")
        boi = st.sidebar.number_input("Inicial Bo", format="%.6f")
        gra = G_method(pr, np,wp, bo, cf, sw0, boi, we,pi,t,salinity)
        fig = px.scatter(gra, x='We*Bw/Et', y='F/Eo+Efw', title='Gráfico Havlena y Odeh')
        fig.update_xaxes(title_text='We*Bw/Et')
        fig.update_yaxes(title_text='F/Eo+Efw')
        x = gra['We*Bw/Et']
        y = gra['F/Eo+Efw']
        slope, intercept, r, p, se = stats.linregress(x, y)
        yt = intercept + (x * slope)
        fig.add_scatter(x=x, y=yt, mode='lines', line=dict(color='green'),
                        name='Regression Lineal')
        text = f"N [MMStb]: {intercept / 1000000:.4f}"
        fig.add_annotation(text=text, x=1e7, y=2e8)
        st.write(text)
        st.plotly_chart(fig)

    except Exception as e:
        st.write("")


elif options== "Analytical method":
    try:
        aq_por = st.sidebar.number_input("Porosity")
        theta = st.sidebar.number_input("Angle")
        k = st.sidebar.number_input("Permeability")
        cf = st.sidebar.number_input("Formation Compressibility", format="%.6f")
        t = st.sidebar.number_input("Temperature")
        salinity = st.sidebar.number_input("Salinity")
        water_visc = st.sidebar.number_input("Water viscosity ")
        se_np = st.selectbox("Production oil cumulate", df.columns)
        se_wp = st.selectbox("Production water cumulate", df.columns)
        se_ppvt = st.selectbox("Pressure PVT", df_pvt.columns)
        se_oil_fvf = st.selectbox("Bo", df_pvt.columns)
        sw0 = st.sidebar.number_input("Saturation water")
        pi = st.sidebar.number_input("Inicial Pressure")
        aq_radius = st.slider("Aquifer Radius", 1000, 20000)
        res_radius = st.slider("Reservoir Radius", 100, 5000)
        aq_thickness = st.slider("Aquifer Thickness", 1, 200)
        N = st.slider("Original Oil in Place",1e+6,200e+6)
        Np = df[se_np].values
        wp = df[se_wp].values
        ppvt_col = se_ppvt
        oil_fvf_col = se_oil_fvf

        P_nueva=calcuted_pressure(Np, wp, cf, t, salinity, df_pvt, aq_radius, res_radius, aq_thickness, aq_por, theta, k, water_visc,
                pi, sw0, N, ppvt_col, oil_fvf_col)

        fig = px.scatter(df, x='Date', y='Pressure', title='Graph P vs t')
        fig.add_scatter(x=df['Date'], y=P_nueva, mode='lines', line=dict(color='green'),
                        name='Presion Calculada')
        fig.update_xaxes(title_text='Tiempo (Años)')
        fig.update_yaxes(title_text='Presion (psia)')
        text = f"N [MMStb]: {N / 1000000:.4f}"
        st.write(text)
        st.plotly_chart(fig)
    except Exception as e:
        st.write("")