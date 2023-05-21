from calendar import monthrange
import datetime
import numpy as np
import pandas as pd


def days_in_month(date):
    if isinstance(date, datetime.datetime):
        return monthrange(date.year, date.month)[1]
    else:
        raise ValueError("Argument is not of type datetime")


def interp_from_dates(date_interp: datetime.datetime, x_dates, y_values,
                      left=None, right=None):
    """
    Interpolation function that accepts dates as x to interpolate between y_values,
    given a new date between x_dates

    Parameters
    ----------
    date_interp: array-like of datetime objects
        the new date to interpolate
    x_dates: array-like of datetime objects
        the x_dates values to use as regression
    y_values: array-like
        the y_values to use as regression
    left: optional float or complex corresponding to y_values
        Value to return for x < x_dates[0], default is y_values[0]
    right: optional float or complex corresponding to y_values
        Value to return for x > x_dates[-1], default is y_values[-1]
    Returns
    -------
    float or ndarray
        an interpolated value or values between y_values

    """

    # Type checking before proceeding with operations
    permitted_arrays = (list, np.ndarray, pd.Series)
    if isinstance(date_interp, permitted_arrays):
        is_array = True
        if not all(isinstance(x, datetime.datetime) for x in date_interp):
            raise ValueError("date_interp should only contain datetime objects")
    else:
        is_array = False
        if not isinstance(date_interp, datetime.datetime):
            raise ValueError(f"{date_interp} is not a datetime object")

    if isinstance(x_dates, permitted_arrays):
        if not all(isinstance(x, datetime.datetime) for x in x_dates):
            raise ValueError("x_dates should only contain datetime objects")
    else:
        raise ValueError("x_dates should be a list, numpy array or pandas Series")

    if not isinstance(y_values, permitted_arrays):
        raise ValueError("y_values should be a list, numpy array or pandas Series")

    # The handling of numeric values in y_values should be handled by np itself in the
    # interpolation function

    # Get the minimum date in x_dates as the reference
    start_date = x_dates.min()
    # Calculate time deltas using x_dates and the start dates, then convert to seconds
    time_deltas = pd.Series(x_dates - start_date).dt.total_seconds()

    if is_array:
        # do the same for date_interp values
        new_time_delta = pd.Series(date_interp - start_date).dt.total_seconds()
    else:
        new_time_delta = (date_interp - start_date).total_seconds()

    return np.interp(new_time_delta, time_deltas, y_values, left=left, right=right)


def interp_dates_row(row, x_result_col, df_input: pd.DataFrame,
                     x_input_col, y_input_col, input_cond_col_name,
                     result_cond_col_name, left=None, right=None):
    """
    A helper function that works on data frame rows using the apply method.
    This function creates an interpolation using interp_from_dates whose input values
    would change based on a condition specified in the target data frame

    Parameters
    ----------
    row: DataFrame
        Usually used with apply method, axis=1
    x_result_col: str
        the x value col name to interpolate in the target row
    df_input: DataFrame
        the input data frame that will be used for regression
    x_input_col: str
        the x column in the input data frame
    y_input_col: str
        the y column in the input data frame
    input_cond_col_name: str
        the condition column in the input data frame
    result_cond_col_name: str
        the condition column in the result data frame
    left: optional float or complex corresponding to y_values
        Value to return for x < x_dates[0], default is y_values[0]
    right: optional float or complex corresponding to y_values
        Value to return for x > x_dates[-1], default is y_values[-1]

    Returns
    -------
    float
        returns an interpolated value of float type
    """

    # Filter the input data frame according to group name
    df_input_filt = df_input.loc[df_input[input_cond_col_name]
                                 == row[result_cond_col_name]]

    if len(df_input_filt) == 0:
        return np.nan
    else:
        return interp_from_dates(row[x_result_col],
                                 df_input_filt[x_input_col],
                                 df_input_filt[y_input_col],
                                 left=left, right=right)


def material_bal_var_type(data, numb_or_column):
    """
    Function to check the data types of the material balance equation terms

    Parameters
        ----------
    data: Pandas Dataframe
        Contains the production information for a single entity
    numb_or_column: Dictionary
        Python dictionary containing the names of some columns of the dataframe in
        order to check their data types

    Returns
    -------
    Pandas Dataframe:
        Returns the original Pandas Dataframe
    """
    # Make a copy of the original dataframe
    df = data.copy()
    if not isinstance(data, pd.DataFrame):
        raise TypeError("The input data should be a pandas dataframe")

    # Define internal names for column in the DataFrame
    for col, arg in numb_or_column.items():
        if isinstance(arg, (int, float)):
            df[col] = arg
            numb_or_column[col] = col
        elif isinstance(arg, str):
            df.rename(columns={arg: col}, inplace=True)
        else:
            raise TypeError(f"{arg} should be either a numeric value or string "
                            f"indicating a column in the DataFrame")

    return df


def material_bal_numerical_data(vector):
    """
    Function to check the numerical data types of the arguments of the material balance
    functions

    Parameters
    ----------
    vector: List or numpy array
        List or array of numerical arguments for each function
        Contains the production information for a single entity

    Returns
    -------
    A message showing if there is presence of TypeError """
    for value in vector:
        if not isinstance(value, (int, float)):
            raise TypeError(f"{value} should be either an int or float")


def variable_type(obj):
    """
        :param obj: variable to be converted as an array, it may be a list or a float
        :return: an array if obj is entered in the right data format
    """
    if isinstance(obj, np.ndarray):
        array = obj
    elif isinstance(obj, list):
        array = np.array(obj)
    elif isinstance(obj, float):
        array = np.array(obj)
    else:
        raise ValueError("Please enter measured values as float, list or array")
    return array
  