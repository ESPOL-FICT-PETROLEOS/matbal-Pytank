import pandas as pd
import pytest


@pytest.fixture()
def prod_data_1():

    data = {"Np": [1000, 2000, 3000],  # bls
            "Gp": [100000, 200000, 300000],  # scf, assuming rs of 100 scf/STB
            "Wp": [100, 200, 300]}  # bls

    df = pd.DataFrame(data)

    yield df


@pytest.fixture()
def press_data_uw_1():

    df = pd.read_csv("data_for_tests/other_data/press_data_1.csv")
    df["Date"] = df["Date"].astype("datetime64")

    yield df
