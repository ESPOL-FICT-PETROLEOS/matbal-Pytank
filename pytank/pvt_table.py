
from pytank import pvt_correlations as pvt
import pandas as pd
import numpy as np


def pvt_table(p_sep, t_sep, api, rsp, sg_sep, tres, den_sto, p_res, salinity, jump,
              units=1) -> pd.DataFrame:

    """The pvt_table data frame is a table that has a size of the length of Pres(it
    is printed regarding a range) vs 6 columns, where within the first column are shown
    the pressure values (from 0 to Pres) whereas within each of the next 6 columns, are
    displayed 6 PVT properties such as the oil density, FVF_oil, GOR, rsw, FVF_water,
    and water compressibility evaluated at the values of P from 0 to the reservoir
    pressure. The required input data is: Separator pressure [psia], Separator
    Temperature [F], Stock tank oil gravity [API], Separator producing GOR [scf/STB],
    Separator gas Specific Gravity, Reservoir Temperature [F], Density of stock-tank oil
    [lb/cu ft], Reservoir Pressure [psia], water salinity, and Jump(it indicates the
    step size of the iteration over the pressure).It also requires the Pb from Velarde
    correlation, which is called within the function.

    Parameters
    ----------
    p_sep: int or float
        Separator pressure [psia]
    t_sep: int or float
        Separator temperature [F]
    api: int or float
        Stock tank oil gravity [API]
    rsp: int or float
        Separator producing GOR [scf/STB]
    sg_sep: int or float
        Separator gas specific gravity
    tres: int or float
        Reservoir temperature [F]
    den_sto: int or float
        Density of stock-tank oil [lb/cu ft]
    p_res: int or float
        Initial reservoir pressure [psia]
    salinity: int or float
        Water Salinity [ppm]
    jump: int or float
        step size of the iteration over the pressure
    units: unit system
        if units = 1, it will use field units. Otherwise, it will be used metric units

    Returns
    -------
    Pandas Dataframe:
        Returns a pandas Dataframe with some pvt properties
    """

    # Call the bubble point pressure function
    pb = pvt.Pb_Velarde(api, sg_sep, tres, p_sep, t_sep, rsp, units)

    # Naming the dataframe columns
    pressure = "pressure[psia]"
    gor = "gor[scf/stb]"
    density = "density[lb/cu ft]"
    oil_fvf = "oil_fvf[rb/stb]"
    rsw = "gas_water_Rs[scf/stb]"
    water_fvf = "water_fvf[rb/stb]"
    water_comp = "water_comp[1/psia]"

    # List containing the names of all columns
    columns = [pressure, gor, density, oil_fvf, rsw, water_fvf, water_comp]

    # Creation of empty dataframes; under, at, and above the bubble pressure point
    df_under = pd.DataFrame(columns=columns)
    df_pb = pd.DataFrame(columns=columns)
    df_above = pd.DataFrame(columns=columns)

    # Loop to iterate from 0 to reservoir pressure
    for p in np.arange(0, p_res+1, jump):

        # Condition when pressure is less than the bubble point pressure
        if (p > 0) and (p < pb):

            gor_velarde2 = pvt.Solution_GOR_Velarde2(sg_sep, api, tres, p, p_sep, t_sep,
                                                     rsp, units)
            den_underpb = pvt.Den_underPb(sg_sep, tres, p, api, p_sep, t_sep, rsp,
                                          units)
            fvf_underpb = pvt.FVF_underPb(den_sto, p_sep, t_sep, api, rsp, sg_sep, tres,
                                          p, units)
            rsw_under = pvt.RS_bw(p, tres, salinity, units)
            water_fvf_under = pvt.Bo_bw(p, tres, salinity, units)
            water_comp_under = pvt.comp_bw_nogas(p, tres, salinity, units)
            p_under = p

            df_under = df_under.append({pressure: p_under, gor: gor_velarde2,
                                            density: den_underpb, oil_fvf: fvf_underpb,
                                            rsw: rsw_under, water_fvf: water_fvf_under,
                                            water_comp: water_comp_under},
                                       ignore_index=True)

        # Condition when pressures is equal to the bubble point pressure
        elif p == pb:
            df_pb = pd.DataFrame({pressure: [p],
                gor:[pvt.Solution_GOR_Pb_ValkoMcCain(p_sep, t_sep, api, rsp, units)],
                density: [pvt.Den_Pb(sg_sep, tres, p_sep, t_sep, api, rsp, units)],
                oil_fvf: [pvt.FVF_Pb(den_sto, p_sep, t_sep, api, rsp, sg_sep, tres,
                                     units)],
                rsw: [pvt.RS_bw(p, tres, salinity, units)], water_fvf:
                [pvt.Bo_bw(p, tres, salinity, units)],
                water_comp: [pvt.comp_bw_nogas(p, tres, salinity, units)]})

        # Condition when pressure is greater than the bubble point pressure
        elif (p > pb) and (p <= p_res):

            den_abovepb = pvt.Den_abovePb(p, sg_sep, tres, p_sep, t_sep, api, rsp,
                                          p_res, units)
            gor_velarde2_above = pvt.Solution_GOR_Pb_ValkoMcCain(p_sep, t_sep, api, rsp,
                                                                 units)
            fvf_abovepb = pvt.FVF_abovePb(p, den_sto, p_sep, t_sep, api, rsp, sg_sep,
                                          tres, p_res, units)
            rsw_above = pvt.RS_bw(p, tres, salinity, units)
            water_fvf_above = pvt.Bo_bw(p, tres, salinity, units)
            water_comp_above = pvt.comp_bw_nogas(p, tres, salinity, units)
            p_above = p

            df_above = df_above.append({pressure: p_above, gor: gor_velarde2_above,
                                            density: den_abovepb, oil_fvf: fvf_abovepb,
                                        rsw: rsw_above, water_fvf: water_fvf_above,
                                        water_comp: water_comp_above},
                                       ignore_index=True)

    # Concatenation of the pvt dataframes regarding their relationship to the bubble
    # point pressure
    pvt_dataframe = df_under.append([df_pb, df_above])

    return pvt_dataframe


