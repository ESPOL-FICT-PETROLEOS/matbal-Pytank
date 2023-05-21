from pytank.material_balance import underground_withdrawal, pressure_vol_avg
import numpy as np
import pytest
import pandas as pd


def test_uw_only_prod(prod_data_1):
    """Test expected calculation only"""
    df = prod_data_1

    oil_fvf = 1.2
    water_fvf = 1
    gas_fvf = 0.001
    rs = 100
    rsw = 0

    uw_calculated = underground_withdrawal(df, "Np", "Wp", "Gp",
                                           oil_fvf, water_fvf, gas_fvf,
                                           rs, rsw)

    uw_expected = np.array([1300, 2600, 3900])

    # Check first if it is a numpy array type
    assert isinstance(uw_calculated, np.ndarray)
    # Check for the result itself
    assert pytest.approx(uw_calculated) == uw_expected


def test_uw_variable_bo(prod_data_1):
    """Test for changing values of a fluid property"""
    df = prod_data_1

    df["Bo"] = [1.2, 1.3, 1.4]
    water_fvf = 1
    gas_fvf = 0.001
    rs = 100
    rsw = 0

    uw_calculated = underground_withdrawal(df, "Np", "Wp", "Gp",
                                           "Bo", water_fvf, gas_fvf,
                                           rs, rsw)

    uw_expected = np.array([1300, 2700, 4200])

    # Check first if it is a numpy array type
    assert isinstance(uw_calculated, np.ndarray)
    # Check for the result itself
    assert pytest.approx(uw_calculated) == uw_expected


def test_uw_negative_gas_values(prod_data_1):
    """Test the case when the gas withdrawal is negative"""
    df = prod_data_1

    oil_fvf = 1.2
    water_fvf = 1
    gas_fvf = 0.001
    rs = 120
    rsw = 0

    with pytest.raises(ArithmeticError):
        _ = underground_withdrawal(df, "Np", "Wp", "Gp",
                                   oil_fvf, water_fvf, gas_fvf,
                                   rs, rsw)


def test_uw_wrong_input_df():
    """Test for the case when the input data is not a DataFrame"""
    df = np.array([1000, 2000, 3000])

    oil_fvf = 1.2
    water_fvf = 1
    gas_fvf = 0.001
    rs = 120
    rsw = 0

    with pytest.raises(TypeError) as exception_info:
        _ = underground_withdrawal(df, "Np", "Wp", "Gp",
                                   oil_fvf, water_fvf, gas_fvf,
                                   rs, rsw)

    exception_expected = "The input data should be a pandas dataframe"

    assert exception_info.match(exception_expected)


def test_uw_wrong_fluid_input(prod_data_1):
    """Test for wrong fluid input value"""
    df = prod_data_1

    oil_fvf = [1.2, 1.3, 1.4]
    water_fvf = 1
    gas_fvf = 0.001
    rs = 120
    rsw = 0

    with pytest.raises(TypeError) as exception_info:
        _ = underground_withdrawal(df, "Np", "Wp", "Gp",
                                   oil_fvf, water_fvf, gas_fvf,
                                   rs, rsw)

    exception_expected = f"{oil_fvf} should be either a numeric value or string " \
                         f"indicating a column in the DataFrame"

    assert exception_info.value.args[0] == exception_expected


def test_uw_wrong_column(prod_data_1):
    """Test for a column that does not exist in the DataFrame"""
    df = prod_data_1

    oil_fvf = 1.2
    water_fvf = 1
    gas_fvf = 0.001
    rs = 120
    rsw = 0

    with pytest.raises(KeyError):
        _ = underground_withdrawal(df, "NpP", "Wp", "Gp",
                                   oil_fvf, water_fvf, gas_fvf,
                                   rs, rsw)


def test_press_avg_data_1(press_data_uw_1):
    """Test basic case using press_data_uw_1 fixture"""
    df = press_data_uw_1
    start_date = df["Date"].min()
    freq = "12MS"

    dates = pd.date_range(start_date, freq=freq, periods=4)
    data = {"Date": dates, "Pressure": [3697.517, 2809.678, np.nan, 1500]}
    result_expected = pd.DataFrame(data)

    result = pressure_vol_avg(df, "Well", "Date", "Pressure", "UW", freq)

    assert (pytest.approx(result["Pressure"], nan_ok=True)
            == result_expected["Pressure"])
    assert all(result_expected["Date"] == result["Date"])


def test_press_avg_data_1_middle(press_data_uw_1):
    """Test basic case using press_data_uw_1 fixture and middle argument"""
    df = press_data_uw_1
    start_date = df["Date"].min()
    freq = "12MS"

    dates = pd.date_range(start_date, freq=freq, periods=4)
    deltas = pd.to_timedelta([366, 365, 365, 365], unit="d") / 2
    new_date = dates + deltas

    data = {"Date": new_date, "Pressure": [3697.517, 2809.678, np.nan, 1500]}
    result_expected = pd.DataFrame(data)

    result = pressure_vol_avg(df, "Well", "Date", "Pressure", "UW", freq, "middle")

    assert (pytest.approx(result["Pressure"], nan_ok=True)
            == result_expected["Pressure"])
    assert all(result_expected["Date"] == result["Date"])


def test_press_avg_data_1_end(press_data_uw_1):
    """Test basic case using press_data_uw_1 fixture and "end" argument"""
    df = press_data_uw_1
    start_date = df["Date"].min()
    freq = "12MS"

    dates = pd.date_range(start_date, freq=freq, periods=4)
    # This is the only difference compared to the middle test
    deltas = pd.to_timedelta([366, 365, 365, 365], unit="d")
    new_date = dates + deltas

    data = {"Date": new_date, "Pressure": [3697.517, 2809.678, np.nan, 1500]}
    result_expected = pd.DataFrame(data)

    result = pressure_vol_avg(df, "Well", "Date", "Pressure", "UW", freq, "end")

    assert (pytest.approx(result["Pressure"], nan_ok=True)
            == result_expected["Pressure"])
    assert all(result_expected["Date"] == result["Date"])
