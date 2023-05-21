# %%
from pytank.aquifer import aquifer_fetkovich
from pytank.aquifer import aquifer_carter_tracy
import pandas as pd
import pytest


# Test on normal arguments
# TODO aquifer fetkovich test for boundary type=constant_pressure and flow_type=radial
# TODO aquifer fetkovich test for boundary type=constant_pressure and flow_type=radial
# TODO boundary type=no_flow and flow_type=linear
# TODO boundary type=constant_pressure and flow_type=linear

# Test aquifer fetkovich function for: boundary type=no_flow and flow_type=radial
# This data was taken from Example 10-10 in the Reservoir Engineering Handbook.
# There you will find the solution.
# Reference: Tarek,A.(2019). Chapter 10. Water Influx,
# Reservoir Engineering Handbook( 5th edition, pages.746-747).
# The are small differences between expected results and calculated results
# due to the number of used decimals that is why this test fail.

def test_on_aquifer_fetkovich_output():
    aq_radius = 46000
    res_radius = 9200
    aq_thickness = 100
    # -fi = 0.25
    # +phi = 0.25
    phi = 0.25
    ct = 0.000007
    pr = [2740, 2500, 2290, 2109, 1949]
    theta = 140
    k = 200
    water_visc = 0.55
    time_step = [0, 365, 730, 1095, 1460]

    df_calculated = aquifer_fetkovich(aq_radius, res_radius, aq_thickness, phi, ct, pr,
                                      theta, k, water_visc, time_step,
                                      boundary_type='no_flow', flow_type='radial',
                                      width=None, length=None)
    df = pd.DataFrame(columns=['Elapsed time', 'Delta We', 'Cumulative We'])
    df['Elapsed time'] = [0, 365, 730, 1095, 1460]
    df['Delta We'] = [0, 3925000, 9615000, 11970000, 12461000]
    df['Cumulative We'] = [0, 3925000, 13540000, 25510000, 37971000]
    df = df.set_index('Elapsed time')
    df_expected = df
    assert isinstance(df_calculated, pd.DataFrame)
    assert pytest.approx(df_calculated) == df_expected


# Test for case td > 100
# Reference: Tarek,A.(2019). Chapter 10. Water Influx,
# Reservoir Engineering Handbook( 5th edition, page. 740)
# Data taken from example 10.6 (page. 693) since example 10.9 does not
# specify the values used for each parameter,
# but the dimensionless time function is the same in both cases.
# However there are differences in water influx constant b which is specified as 20.4
# in example 10.9 but in example 10.6, with the given values
# the calculated b is 22.4 that's why this test fail.
# Also we noticed there is an error in calculations shown in example 10.9.
def test_on_aquifer_carter_tracy_output():
    aq_por = 0.2
    ct = 0.000001
    res_radius = 2000
    aq_thickness = 25
    theta = 360
    aq_perm = 100
    water_visc = 0.8
    pr = [2500, 2490, 2472, 2444, 2408]
    time = [0, 182.5, 365.0, 547.5, 730.0]

    df_calculated = aquifer_carter_tracy(aq_por, ct, res_radius, aq_thickness, theta,
                                         aq_perm,
                                         water_visc, pr, time)
    df = pd.DataFrame(columns=['Elapsed time, days', 'Cumulative water influx, bbl'])
    df['Elapsed time, days'] = [0, 182.5, 365.0, 547.5, 730.0]
    df['Cumulative water influx, bbl'] = [0, 12265.82, 44545.75, 106302.54, 204312.42]
    df = df.set_index('Elapsed time, days')
    df_expected = df
    assert isinstance(df_calculated, pd.DataFrame)
    assert pytest.approx(df_calculated) == df_expected


# TODO Test for case td < 100
# There is an example to test for this case in the L.P Dake book,
# The practice of reservoir engineering but unfortunately the data
# to be used is not specified and also, the approach to estimate
# the dimensionless pressure and the dimensionless pressure derivative
# is different from that specified in the Reservoir engineer handbook from Tarek, A.

# Test on bad arguments


def test_on_negative_pressure_fetkovich():
    aq_radius = 46000
    res_radius = 9200
    aq_thickness = 100
    # -fi = 0.25
    # +phi = 0.25
    aq_por = 0.25
    ct = 0.000007
    pr = [2740, 2500, 2290, 2109, -1949]
    theta = 140
    k = 200
    water_visc = 0.55
    time_step = [0, 365, 730, 1095, 1460]
    with pytest.raises(ValueError) as exception_info:
        aquifer_fetkovich(aq_radius, res_radius, aq_thickness, aq_por, ct, pr,
                          theta, k, water_visc, time_step,
                          boundary_type='no_flow', flow_type='radial',
                          width=None, length=None)
    assert exception_info.match("Pressure must be greater than zero")


def test_on_negative_pressure_carter_tracy():
    aq_por = 0.2
    ct = 0.000001
    res_radius = 2000
    aq_thickness = 25
    theta = 360
    aq_perm = 100
    water_visc = 0.8
    pr = [2500, 2490, 2472, 2444, -2408]
    time = [0, 182.5, 365.0, 547.5, 730.0]
    with pytest.raises(ValueError) as exception_info:
        aquifer_carter_tracy(aq_por, ct, res_radius, aq_thickness, theta, aq_perm,
                             water_visc, pr, time)
    assert exception_info.match("Pressure must be greater than zero")


def test_on_pressure_descendent_order_fetkovich():
    aq_radius = 46000
    res_radius = 9200
    aq_thickness = 100
    # -fi = 0.25
    # +phi = 0.25
    aq_por = 0.25
    ct = 0.000007
    pr = [1949, 2109, 2290, 2500, 2740]
    theta = 140
    k = 200
    water_visc = 0.55
    time_step = [0, 365, 730, 1095, 1460]
    with pytest.raises(ValueError) as exception_info:
        aquifer_fetkovich(aq_radius, res_radius, aq_thickness, aq_por, ct, pr,
                          theta, k, water_visc, time_step,
                          boundary_type='no_flow', flow_type='radial',
                          width=None, length=None)
    assert exception_info.match("Pressure array must be in descendant order")


def test_on_pressure_descendent_order_carter_tracy():
    aq_por = 0.2
    ct = 0.000001
    res_radius = 2000
    aq_thickness = 25
    theta = 360
    aq_perm = 100
    water_visc = 0.8
    pr = [2408, 2444, 2472, 2490, 2500]
    time = [0, 182.5, 365.0, 547.5, 730.0]
    with pytest.raises(ValueError) as exception_info:
        aquifer_carter_tracy(aq_por, ct, res_radius, aq_thickness, theta, aq_perm,
                             water_visc, pr, time)
    assert exception_info.match("Pressure array must be in descendant order")


def test_on_pressure_time_array_dimensions_fetkovich():
    aq_radius = 46000
    res_radius = 9200
    aq_thickness = 100
    # -fi = 0.25
    # +phi = 0.25
    aq_por = 0.25
    ct = 0.000007
    pr = [2740, 2500, 2290, 2109, 1949]
    theta = 140
    k = 200
    water_visc = 0.55
    time_step = [0, 365, 730, 1095]
    with pytest.raises(ValueError) as exception_info:
        aquifer_fetkovich(aq_radius, res_radius, aq_thickness, aq_por, ct, pr,
                          theta, k, water_visc, time_step,
                          boundary_type='no_flow', flow_type='radial',
                          width=None, length=None)
    assert exception_info.match("Dimensions of pressure array and time array "
                                "should be equal,"
                                " please verify your input")


def test_on_pressure_time_array_dimensions_carter_tracy():
    aq_por = 0.2
    ct = 0.000001
    res_radius = 2000
    aq_thickness = 25
    theta = 360
    aq_perm = 100
    water_visc = 0.8
    pr = [2500, 2490, 2472, 2444, 2408]
    time = [0, 182.5, 365.0, 547.5]
    with pytest.raises(ValueError) as exception_info:
        aquifer_carter_tracy(aq_por, ct, res_radius, aq_thickness, theta, aq_perm,
                             water_visc, pr, time)
    assert exception_info.match("Dimensions of pressure array and time array "
                                "should be equal,"
                                " please verify your input")


def test_on_pressure_type_input_fetkovich():
    aq_radius = 46000
    res_radius = 9200
    aq_thickness = 100
    # -fi = 0.25
    # +phi = 0.25
    aq_por = 0.25
    ct = 0.000007
    pr = "2740, 2500, 2290, 2109, 1949"
    theta = 140
    k = 200
    water_visc = 0.55
    time_step = [0, 365, 730, 1095, 1460]
    with pytest.raises(ValueError) as exception_info:
        aquifer_fetkovich(aq_radius, res_radius, aq_thickness, aq_por, ct, pr,
                          theta, k, water_visc, time_step,
                          boundary_type='no_flow', flow_type='radial',
                          width=None, length=None)
    assert exception_info.match("Please enter measured values as float, list or array")


def test_on_pressure_type_input_carter_tracy():
    aq_por = 0.2
    ct = 0.000001
    res_radius = 2000
    aq_thickness = 25
    theta = 360
    aq_perm = 100
    water_visc = 0.8
    pr = "2500, 2490, 2472, 2444, 2408"
    time = [0, 182.5, 365.0, 547.5]
    with pytest.raises(ValueError) as exception_info:
        aquifer_carter_tracy(aq_por, ct, res_radius, aq_thickness, theta, aq_perm,
                             water_visc, pr, time)
    assert exception_info.match("Please enter measured values as float, list or array")


def test_on_time_type_input_fetkovich():
    aq_radius = 46000
    res_radius = 9200
    aq_thickness = 100
    # -fi = 0.25
    # +phi = 0.25
    aq_por = 0.25
    ct = 0.000007
    pr = [2740, 2500, 2290, 2109, 1949]
    theta = 140
    k = 200
    water_visc = 0.55
    time_step = "360"
    with pytest.raises(ValueError) as exception_info:
        aquifer_fetkovich(aq_radius, res_radius, aq_thickness, aq_por, ct, pr,
                          theta, k, water_visc, time_step,
                          boundary_type='no_flow', flow_type='radial',
                          width=None, length=None)
    assert exception_info.match("Please enter measured values as float, list or array")


def test_on_time_type_input_carter_tracy():
    aq_por = 0.2
    ct = 0.000001
    res_radius = 2000
    aq_thickness = 25
    theta = 360
    aq_perm = 100
    water_visc = 0.8
    pr = [2500, 2490, 2472, 2444, 2408]
    time = "512"
    with pytest.raises(ValueError) as exception_info:
        aquifer_carter_tracy(aq_por, ct, res_radius, aq_thickness, theta, aq_perm,
                             water_visc, pr, time)
    assert exception_info.match("Please enter measured values as float, list or array")


def test_on_flow_type_input():
    aq_radius = 46000
    res_radius = 9200
    aq_thickness = 100
    # -fi = 0.25
    # +phi = 0.25
    aq_por = 0.25
    ct = 0.000007
    pr = [2740, 2500, 2290, 2109, 1949]
    theta = 140
    k = 200
    water_visc = 0.55
    time_step = [0, 365, 730, 1095, 1460]
    with pytest.raises(ValueError) as exception_info:
        aquifer_fetkovich(aq_radius, res_radius, aq_thickness, aq_por, ct, pr,
                          theta, k, water_visc, time_step,
                          boundary_type='no_flow', flow_type='linear',
                          width=None, length=None)
    assert exception_info.match("When using linear flow, "
                                "width and length are required arguments")
