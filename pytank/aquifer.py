# File to code aquifer-related functions or classes
from pytank.utilities import variable_type
import math
import numpy as np
import pandas as pd


def schilthuis_aq(pi, p, k, h, uw, ra, re):
    """
    This function calculates the aquifer water influx according to schilthuis equation
    :param pi: Initial reservoir pressure (psi)
    :param p: Current reservoir pressure, might be scalar, list or numpy array (psi)
    :param k: effective permeability to water (mili darcy)
    :param h: Estimated reservoir thickness in the aquifer zone (ft)
    :param uw: Water viscosity (cp)
    :param ra: aquifer radius (ft)
    :param re: hydrocarbon zone radius (ft)
    :return: a single or array-like value of aquifer influx (bbl)
    """

    # Convert to numpy array in case the user has input a scalar or list
    # TODO Check only for numerical values and lists, otherwise throw an Exception
    if isinstance(p, np.ndarray):
        p_np = p
    else:
        p_np = np.array(p)

    # Calculate water influx (Field units)
    # TODO handle different units
    c = 0.00708 * k * h / (uw * np.log(ra / re))
    ew = c * (pi - p_np)

    return ew


def aquifer_fetkovich(aq_radius, res_radius, aq_thickness, aq_por, ct, pr, theta, k,
                      water_visc, time_step,
                      boundary_type='no_flow', flow_type='radial',
                      width=None, length=None):
    """
        To estimate water influx using the Fetkovich's method we need to estimate:
        Wi, Wei and J.
        :param aq_radius: radius of the aquifer, ft
        :param res_radius: radius of the reservoir, ft
        :param aq_thickness: thickness of the aquifer, ft
        #param fi: changed by phi
        :param aq_por: porosity of the aquifer
        :param ct: total compressibility coefficient, psi-1
        :param pr: measured reservoir pressure, may be an integer, a float, list or
        numpy array, psi
        :param theta: encroachment angle
        :param k: permeability of the aquifer, md
        :param water_visc: viscosity of water, cp
        :param time_step: time step may be an integer, a float, list or
        numpy array, days
        :param boundary_type: default value = 'no_flow',
        options are: 'no_flow', 'constant_pressure', 'infinite'
        :param flow_type: default value = 'radial', options are: 'radial', 'linear'
        :param width: width of the linear aquifer,
        * parameter required only for linear flow, ft
        :param length: length of the linear aquifer,
        * parameter required only for linear flow, ft
        :return: a DataFrame containing the cumulative water influx, bbl
    """
    if flow_type == 'linear' and (width is None or length is None):
        raise ValueError("When using linear flow, "
                         "width and length are required arguments")
    # Check if pressure and time step are arrays, list or floats
    pr_array = variable_type(pr)
    delta_t = variable_type(time_step)
    # Check if pressure array is not in descendant order throw an error
    order = pd.Series(pr_array).is_monotonic_decreasing
    if order is False:
        raise ValueError("Pressure array must be in descendant order")
    # Check if time step and pressure dimensions are equal
    # this can be done if time step is entered as array
    if not all(pr_array > 0):
        raise ValueError("Pressure must be greater than zero")
    if isinstance(pr_array, np.ndarray) and isinstance(delta_t, np.ndarray):
        dim_pr = np.size(pr_array)
        dim_time = np.size(delta_t)
        if dim_pr != dim_time:
            raise ValueError("Dimensions of pressure array and time array "
                             "should be equal,"
                             " please verify your input")
    # Calculate the initial volume of water in the aquifer (Wi)
    wi = (math.pi / 5.615) * (aq_radius ** 2 - res_radius ** 2) * aq_thickness * aq_por
    # Calculate the maximum possible water influx (Wei)
    f = theta / 360
    wei = ct * wi * pr_array[0] * f
    # Calculate the aquifer productivity index
    # based on the boundary_type conditions and aquifer geometry (J)
    rd = aq_radius / res_radius

    if boundary_type == "no_flow" and flow_type == "radial":
        j = (0.00708 * k * aq_thickness * f) / (water_visc * (math.log(rd) - 0.75))
    elif boundary_type == "constant_pressure" and flow_type == "radial":
        j = (0.00708 * k * aq_thickness * f) / (water_visc * math.log(rd))
    elif boundary_type == "no_flow" and flow_type == "linear":
        j = (0.003381 * k * width * aq_thickness) / (water_visc * length)
    elif boundary_type == "constant_pressure" and flow_type == "linear":
        j = (0.001127 * k * width * aq_thickness) / (water_visc * length)
    # Calculate the incremental water influx (We)n during the nth time interval
    # Calculate cumulative water influx
    cum_water_influx = 0
    pr = pr_array[0]
    # Average aquifer pressure after removing We bbl of water from the aquifer
    pa = pr_array[0]
    elapsed_time = np.empty((1, 0))
    time_steps = np.array(0)
    df = pd.DataFrame(columns=['Delta We'])
    for ip in range(len(pr_array)):
        pr_avg = (pr + pr_array[ip]) / 2
        if isinstance(delta_t, np.ndarray):
            diff_pr = np.diff(delta_t)
            time_steps = np.append(time_steps, diff_pr)
            we = (wei / pr_array[0]) * \
                 (1 - math.exp((-1 * j * pr_array[0] * time_steps[ip]) / wei)) * \
                 (pa - pr_avg)
            elapsed_time = delta_t
        else:
            we = (wei / pr_array[0]) * \
                 (1 - math.exp((-1 * j * pr_array[0] * delta_t) / wei)) * \
                 (pa - pr_avg)
            elapsed_time = np.append(elapsed_time, delta_t * ip)
        pr = pr_array[ip]
        cum_water_influx = cum_water_influx + we
        pa = pr_array[0] * (1 - (cum_water_influx / wei))
        df = df.append(
            {'Delta We': we, 'Cumulative We': cum_water_influx},
            ignore_index=True)
    df['Elapsed time'] = elapsed_time
    df = df.set_index('Elapsed time')
    return df


def aquifer_carter_tracy(aq_por, ct, res_radius, aq_thickness, theta, aq_perm,
                         water_visc, pr, time):
    """

    :param aq_por: porosity of the aquifer
    :param ct: total compressibility coefficient, psi-1
    :param res_radius: radius of the reservoir, ft
    :param aq_thickness: thickness of the aquifer, ft
    :param theta: encroachment angle
    :param aq_perm: permeability of the aquifer, md
    :param water_visc: viscosity of water, cp
    :param pr: measured reservoir pressure, may be an integer, a float,
    list or numpy array, psi
    :param time: time, may be an integer, a float, list or
    numpy array, days
    :return: a DataFrame containing the cumulative water influx, bbl
    """
    # Check if pressure and time are arrays, lists or floats
    pr_array = variable_type(pr)
    t_array = variable_type(time)
    # Check if pressure array is not in descendant order throw an error
    order = pd.Series(pr_array).is_monotonic_decreasing
    if order is False:
        raise ValueError("Pressure array must be in descendant order")
    if not all(pr_array > 0):
        raise ValueError("Pressure must be greater than zero")
    # Check if time step and pressure dimensions are equal
    # this can be done if time step is entered as array
    if isinstance(pr_array, np.ndarray) and isinstance(t_array, np.ndarray):
        dim_pr = np.size(pr_array)
        dim_time = np.size(t_array)
        if dim_pr != dim_time:
            raise ValueError("Dimensions of pressure array and time array "
                             "should be equal,"
                             " please verify your input")
    # Calculate the van Everdingen-Hurst water influx constant
    f = theta / 360
    b = 1.119 * aq_por * ct * (res_radius ** 2) * aq_thickness * f

    # Estimate dimensionless time (tD)
    cte = 0.006328 * aq_perm / (aq_por * water_visc * ct * (res_radius ** 2))
    td = np.where(t_array > 0, t_array * cte, 0)

    # Calculate the total pressure drop (Pi-Pn) as an array, for each time step n.
    pr_drop = np.where(pr_array > 0, pr_array[0] - pr_array, 1)

    # Estimate the dimensionless pressure
    pr_d = np.where(td > 100, 0.5 * (np.log(td) + 0.80907),
                    ((370.529 * np.sqrt(td)) + (137.582 * td) + (5.69549 * (
                            td ** 1.5))) / (
                            328.834 + (265.488 * np.sqrt(td)) + (45.2157 * td) +
                            (td ** 1.5)))
    # Estimate the dimensionless pressure derivative
    e = 716.441 + (46.7984 * (td ** 0.5)) + (270.038 * td) + (71.0098 * (td ** 1.5))
    d = (1296.86 * (td ** 0.5)) + (1204.73 * td) + \
        (618.618 * (td ** 1.5)) + (538.072 * (td ** 2)) + \
        (142.41 * (td ** 2.5))
    pr_deriv = np.where(td > 100, 1 / (2 * td), e / d)
    print(pr_deriv)

    # Calculate the cumulative water influx at any time, ti
    df = pd.DataFrame(columns=['Cumulative water influx, bbl'])
    we = 0
    df = df.append(
            {'Cumulative water influx, bbl': we},
            ignore_index=True)
    for i in np.arange(1, len(td)):

        a1 = td[i] - td[i-1]
        a2 = b * pr_drop[i]
        a3 = we * pr_deriv[i]
        a4 = pr_d[i]
        a5 = td[i-1] * pr_deriv[i]
        cum_influx_water = we + (a1 * ((a2 - a3) / (a4 - a5)))
        we = cum_influx_water
        df = df.append(
            {'Cumulative water influx, bbl': we},
            ignore_index=True)
    df['Elapsed time, days'] = t_array
    df = df.set_index('Elapsed time, days')
    return df
