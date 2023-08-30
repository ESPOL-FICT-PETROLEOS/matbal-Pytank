# File to code pvt-related functions or classes
import pandas as pd
from scipy import interpolate
def interp_pvt_matbal(pvt:pd.DataFrame, press_col_name, prop_col_name, press_target):
    """
    :param pvt: DataFrame for pvt data
    :param press_col_name: String indicating the column name of pressure values
    :param prop_col_name: String indicating the column name of pvt property values to be interpolated
    :param press_target: Numeric value indicating the reservoir pressure target value to interpolate to
    :return: Interpolated pvt property value that match with the pressure target
    """
    x = pvt[press_col_name]
    y = pvt[prop_col_name]
    function = interpolate.interp1d(x, y, fill_value='extrapolate')
    return function(press_target)
