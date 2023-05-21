import numpy as np

"""----------------------------Oil Correlations--------------------------------------"""

def Solution_GOR_Pb_ValkoMcCain(Psep, Tsep, API, Rsp, Units=1):
    """The ValkoMcCain correlation calculates the solution GOR at bubblepoint pressure,
    the required input data is: Separator pressure [psia], Separator Temperature [F],
    Separator producing GOR [scf/STB] and Stock tank oil gravity [API]

    If Units=1 use:
        Pressure in psia
        Temperature in F
        Stock-Tank oil gravity in API
        Solution Gas Oil Ratio @ Pb in Scf/STB

    If Units=2 use:
        Pressure in MPa
        Temperature in C
        Stock-Tank oil gravity in API
        Solution Gas Oil Ratio @ Pb in m3/m3"""

    if Units == 1:

        Psep_array = np.array(Psep)
        Tsep_array = np.array(Tsep)
        API_array = np.array(API)
        Rsp_array = np.array(Rsp)

    else:
        Psep_array = np.array(Psep)
        Psep_array = Psep_array / (1.45037738 * 10 ** 2)
        Tsep_array = np.array(Tsep)
        Tsep_array = (Tsep_array * 9 / 5) + 32
        API_array = np.array(API)
        Rsp_array = np.array(Rsp)
        Rsp_array = Rsp_array * (5.61458333)

    if Psep_array > 0 and Tsep_array > 0 and API_array > 0 and Rsp_array > 0:

        lnPsep = np.log(Psep_array)
        lnTsep = np.log(Tsep_array)

        #              VARn                C0n            C1n             C2n
        Zntable = {"VAR1": lnPsep, "C01": -8.005, "C11": 2.7, "C21": -0.161,
                   "VAR2": lnTsep, "C02": 1.224, "C12": -0.5, "C22": 0,
                   "VAR3": API_array, "C03": -1.587, "C13": 0.0441, "C23": -2.29e-5}

        z1 = Zntable["C01"] + (Zntable["C11"] * Zntable["VAR1"]) + (
                    Zntable["C21"] * Zntable["VAR1"] ** 2)
        z2 = Zntable["C02"] + (Zntable["C12"] * Zntable["VAR2"]) + (
                    Zntable["C22"] * Zntable["VAR2"] ** 2)
        z3 = Zntable["C03"] + (Zntable["C13"] * Zntable["VAR3"]) + (
                    Zntable["C23"] * Zntable["VAR3"] ** 2)

        sumZ = z1 + z2 + z3
        Rst = np.exp(3.955 + 0.83 * sumZ - 0.024 * (sumZ ** 2) + (0.075 * (sumZ ** 3)))
        Rsb = Rst + Rsp

    elif Psep_array == 0 or Tsep_array == 0 or API_array == 0 and Rsp_array > 0:

        Rsb = 1.1618 * Rsp_array

    else:
        Rsb = "Incorrect Data"

    if Units == 1:

        Rsb = Rsb

    else:
        Rsb = Rsb / (1.78107607 * 10 ** -1)

    return Rsb


# print("Rsb =",Solution_GOR_Pb_ValkoMcCain (Psep, Tsep, API, Rsp))

def Spec_grav_st_ValkoMcCain2(Psep, Tsep, API, Rsp, SGsep, Units=1):
    """ The Spec_grav_st_ValkoMcCain2 correlation calculates the weighted-average
    Specific Gravities of surface gases at Stock-Tank conditions, the required imput
    data is: Separator pressure [psia], Separator Temperature [F], Separator producing
    GOR [scf/STB], Separator gas Specific Gravity and Stock tank oil gravity [API]

        If Units=1 use:
            Pressure in psia
            Temperature in F
            Stock-Tank oil gravity in API
            Solution Gas Oil Ratio @ Pb in Scf/STB

        If Units=2 use:
            Pressure in MPa
            Temperature in C
            Stock-Tank oil gravity in API
            Solution Gas Oil Ratio @ Pb in m3/m3"""

    if Units == 1:

        Psep_array = np.array(Psep)
        Tsep_array = np.array(Tsep)
        API_array = np.array(API)
        Rsp_array = np.array(Rsp)
        SGsep_array = np.array(SGsep)

    else:
        Psep_array = np.array(Psep)
        Psep_array = Psep_array / (1.45037738 * 10 ** 2)
        Tsep_array = np.array(Tsep)
        Tsep_array = (Tsep_array * 9 / 5) + 32
        API_array = np.array(API)
        Rsp_array = np.array(Rsp)
        Rsp_array = Rsp_array * (5.61458333)
        SGsep_array = np.array(SGsep)

    if Psep_array > 0 and Tsep_array > 0 and API_array > 0 and Rsp_array > 0 and \
            SGsep_array > 0:

        lnPsep = np.log(Psep_array)
        lnRsp = np.log(Rsp_array)

        #           VARn            C0n            C1n             C2n    C3n               C4n
        Zntable = {"VAR1": lnPsep, "C01": -17.275, "C11": 7.9597, "C21": -1.1013,
                   "C31": 2.7735e-2, "C41": 3.2287e-3,
                   "VAR2": lnRsp, "C02": -0.3354, "C12": -0.3346, "C22": 0.1956,
                   "C32": -3.4374e-2, "C42": 2.08e-3,
                   "VAR3": API_array, "C03": 3.705, "C13": -0.4273, "C23": 1.818e-2,
                   "C33": -3.459e-4, "C43": 2.505e-6,
                   "VAR4": SGsep_array, "C04": -155.52, "C14": 629.61, "C24": -957.38,
                   "C34": 647.57, "C44": -163.26,
                   "VAR5": Tsep_array, "C05": 2.085, "C15": -7.097e-2, "C25": 9.859e-4,
                   "C35": -6.312e-6, "C45": 1.4e-8}

        z1 = Zntable["C01"] + (Zntable["C11"] * Zntable["VAR1"]) + (
                    Zntable["C21"] * Zntable["VAR1"] ** 2) + (
                     Zntable["C31"] * Zntable["VAR1"] ** 3) + (
                         Zntable["C41"] * Zntable["VAR1"] ** 4)
        z2 = Zntable["C02"] + (Zntable["C12"] * Zntable["VAR2"]) + (
                    Zntable["C22"] * Zntable["VAR2"] ** 2) + (
                     Zntable["C32"] * Zntable["VAR2"] ** 3) + (
                         Zntable["C42"] * Zntable["VAR2"] ** 4)
        z3 = Zntable["C03"] + (Zntable["C13"] * Zntable["VAR3"]) + (
                    Zntable["C23"] * Zntable["VAR3"] ** 2) + (
                     Zntable["C33"] * Zntable["VAR3"] ** 3) + (
                         Zntable["C43"] * Zntable["VAR3"] ** 4)
        z4 = Zntable["C04"] + (Zntable["C14"] * Zntable["VAR4"]) + (
                    Zntable["C24"] * Zntable["VAR4"] ** 2) + (
                     Zntable["C34"] * Zntable["VAR4"] ** 3) + (
                         Zntable["C44"] * Zntable["VAR4"] ** 4)
        z5 = Zntable["C05"] + (Zntable["C15"] * Zntable["VAR5"]) + (
                    Zntable["C25"] * Zntable["VAR5"] ** 2) + (
                     Zntable["C35"] * Zntable["VAR5"] ** 3) + (
                         Zntable["C45"] * Zntable["VAR5"] ** 4)
        sumZ = z1 + z2 + z3 + z4 + z5

        SGst = 1.219 + 0.198 * sumZ + 0.0845 * (sumZ ** 2) + (0.03 * (sumZ ** 3)) + (
                    0.003 * (sumZ ** 4))

    else:
        SGst = "Incorrect Data"

    return SGst


# print("SGst =",Spec_grav_st_ValkoMcCain2(Psep, Tsep, API, Rsp, SGsep))

def Spec_grav_ValkoMcCain2(Psep, Tsep, API, Rsp, SGsep, Units=1):
    """ The ValkoMcCain2 correlation calculates the weighted-average Specific Gravities
    of surface gases, the required imput data is: Separator pressure [psia], Separator
    Temperature [F], Separator producing GOR [scf/STB], Separator gas Specific Gravity
    and Stock tank oil gravity [API]

        If Units=1 use:
            Pressure in psia
            Temperature in F
            Stock-Tank oil gravity in API
            Solution Gas Oil Ratio @ Pb in Scf/STB

        If Units=2 use:
            Pressure in MPa
            Temperature in C
            Stock-Tank oil gravity in API
            Solution Gas Oil Ratio @ Pb in m3/m3"""

    if Units == 1:

        Psep_array = np.array(Psep)
        Tsep_array = np.array(Tsep)
        API_array = np.array(API)
        Rsp_array = np.array(Rsp)
        SGsep_array = np.array(SGsep)

    else:

        Psep_array = np.array(Psep)
        Psep_array = Psep_array / (1.45037738 * 10 ** 2)
        Tsep_array = np.array(Tsep)
        Tsep_array = (Tsep_array * 9 / 5) + 32
        API_array = np.array(API)
        Rsp_array = np.array(Rsp)
        Rsp_array = Rsp_array * (5.61458333)
        SGsep_array = np.array(SGsep)

    if Psep_array > 0 and Tsep_array > 0 and API_array > 0 and Rsp_array > 0 and \
            SGsep_array > 0:

        lnPsep = np.log(Psep_array)
        lnTsep = np.log(Tsep_array)

        #              VARn                C0n            C1n             C2n
        Zntable = {"VAR1": lnPsep, "C01": -8.005, "C11": 2.7, "C21": -0.161,
                   "VAR2": lnTsep, "C02": 1.224, "C12": -0.5, "C22": 0,
                   "VAR3": API_array, "C03": -1.587, "C13": 0.0441, "C23": -2.29e-5}

        z1 = Zntable["C01"] + (Zntable["C11"] * Zntable["VAR1"]) + (
                    Zntable["C21"] * Zntable["VAR1"] ** 2)
        z2 = Zntable["C02"] + (Zntable["C12"] * Zntable["VAR2"]) + (
                    Zntable["C22"] * Zntable["VAR2"] ** 2)
        z3 = Zntable["C03"] + (Zntable["C13"] * Zntable["VAR3"]) + (
                    Zntable["C23"] * Zntable["VAR3"] ** 2)

        sumZ = z1 + z2 + z3
        Rst = np.exp(3.955 + 0.83 * sumZ - 0.024 * (sumZ ** 2) + (0.075 * (sumZ ** 3)))

        SG = ((SGsep_array * Rsp_array) + (
                    Spec_grav_st_ValkoMcCain2(Psep, Tsep, API, Rsp, SGsep) * Rst)) / (
                         Rsp_array + Rst)

    elif Psep_array == 0 or Tsep_array == 0 or API_array == 0 or Rsp_array == 0 and \
            SGsep_array > 0:

        SG = 1.066 * SGsep_array

    else:
        SG = "Incorrect Data"

    return SG


# print("SG =",Spec_grav_ValkoMcCain2 (Psep, Tsep, API, Rsp, SGsep))

def Pb_Velarde(API, SGsep, Tres, Psep, Tsep, Rsp, Units=1):
    """The Velarde correlation calculates the Bubblepoint Pressures at Reservoir
    Temperatures, the required imput data is: Reservoir Temperature [F], Separator
    gas Specific Gravity and Stock tank oil gravity [API]. Also requires the Rsb from
    the ValkoMcCain correlation

        If Units=1 use:
            Pressure in psia
            Temperature in F
            Stock-Tank oil gravity in API
            Solution Gas Oil Ratio @ Pb in Scf/STB

        If Units=2 use:
            Pressure in MPa
            Temperature in C
            Stock-Tank oil gravity in API
            Solution Gas Oil Ratio @ Pb in m3/m3"""

    if Units == 1:

        API_array = np.array(API)
        SGsep_array = np.array(SGsep)
        Tres_array = np.array(Tres)

    else:

        Tres_array = np.array(Tres)
        Tres_array = (Tres_array * 9 / 5) + 32
        API_array = np.array(API)
        SGsep_array = np.array(SGsep)

    if API_array > 0 and SGsep_array > 0 and Tres_array > 0:

        lnRsb = np.log(Solution_GOR_Pb_ValkoMcCain(Psep, Tsep, API, Rsp))

        #              VARn                   C0n              C1n                C2n               C3n
        Zntable = {"VAR1": lnRsb, "C01": -5.48, "C11": -0.0378, "C21": 0.281,
                   "C31": -0.0206,
                   "VAR2": API_array, "C02": 1.27, "C12": -0.0449, "C22": 4.36e-4,
                   "C32": -4.76e-6,
                   "VAR3": SGsep_array, "C03": 4.51, "C13": -10.84, "C23": 8.39,
                   "C33": -2.34,
                   "VAR4": Tres_array, "C04": -0.7835, "C14": 6.23e-3, "C24": -1.22e-5,
                   "C34": 1.03e-8}

        z1 = Zntable["C01"] + (Zntable["C11"] * Zntable["VAR1"]) + (
                    Zntable["C21"] * Zntable["VAR1"] ** 2) + (
                     Zntable["C31"] * Zntable["VAR1"] ** 3)
        z2 = Zntable["C02"] + (Zntable["C12"] * Zntable["VAR2"]) + (
                    Zntable["C22"] * Zntable["VAR2"] ** 2) + (
                     Zntable["C32"] * Zntable["VAR2"] ** 3)
        z3 = Zntable["C03"] + (Zntable["C13"] * Zntable["VAR3"]) + (
                    Zntable["C23"] * Zntable["VAR3"] ** 2) + (
                     Zntable["C33"] * Zntable["VAR3"] ** 3)
        z4 = Zntable["C04"] + (Zntable["C14"] * Zntable["VAR4"]) + (
                    Zntable["C24"] * Zntable["VAR4"] ** 2) + (
                     Zntable["C34"] * Zntable["VAR4"] ** 3)
        sumZ = z1 + z2 + z3 + z4

        Pb = np.exp(7.475 + 0.713 * sumZ + 0.0075 * (sumZ ** 2))

    else:
        Pb = "Missing Data"

    if Units == 1:
        Pb = Pb

    else:
        Pb = Pb * (6.89475729 * 10 ** -3)

    return Pb


# print ("Pb =",Pb_Velarde (API, SGsep, Tres, Psep, Tsep, Rsp))

def Solution_GOR_Velarde2(SGsep, API, Tres, P, Psep, Tsep, Rsp, Units=1):
    """The Velarde2 correlation calculates the Solution Gas-Oil ratios at reservoir
    reservoir pressures less than bubblepoint pressure, the required imput data is:
    Separator gas Specific Gravity, Stock tank oil gravity [API], Reservoir Temperature
    [F] and Pressure [psia]. Also requires the Pb from the Velarde correlation and Rsb
    from the ValkoMcCain correlation

        If Units=1 use:
            Pressure in psia
            Temperature in F
            Stock-Tank oil gravity in API
            Solution Gas Oil Ratio @ Pb in Scf/STB

        If Units=2 use:
            Pressure in MPa
            Temperature in C
            Stock-Tank oil gravity in API
            Solution Gas Oil Ratio @ Pb in m3/m3"""

    if Units == 1:

        SGsep_array = np.array(SGsep)
        API_array = np.array(API)
        Tres_array = np.array(Tres)
        P_array = np.array(P)

    else:
        Tres_array = np.array(Tres)
        Tres_array = (Tres_array * 9 / 5) + 32
        API_array = np.array(API)
        SGsep_array = np.array(SGsep)
        P_array = np.array(P)
        P_array = P_array / (1.45037738 * 10 ** 2)

    if API_array > 0 and SGsep_array > 0 and Tres_array > 0 and P_array < Pb_Velarde(
            API, SGsep, Tres, Psep, Tsep, Rsp):

        Pr = (P_array - 14.7) / (Pb_Velarde(API, SGsep, Tres, Psep, Tsep, Rsp) - 14.7)

        #              An                Bn               Cn
        antable = {"A0": 9.73e-7, "B0": 0.022339, "C0": 0.725167,
                   "A1": 1.672608, "B1": -1.00475, "C1": -1.485480,
                   "A2": 0.92987, "B2": 0.337711, "C2": -0.164741,
                   "A3": 0.247235, "B3": 0.132795, "C3": -0.09133,
                   "A4": 1.056052, "B4": 0.302065, "C4": 0.047094}

        a1 = antable["A0"] * SGsep_array ** antable["A1"] * API_array ** antable[
            "A2"] * Tres_array ** antable["A3"] * (
                     Pb_Velarde(API, SGsep, Tres, Psep, Tsep, Rsp) - 14.7) ** antable[
                 "A4"]
        a2 = antable["B0"] * SGsep_array ** antable["B1"] * API_array ** antable[
            "B2"] * Tres_array ** antable["B3"] * (
                     Pb_Velarde(API, SGsep, Tres, Psep, Tsep, Rsp) - 14.7) ** antable[
                 "B4"]
        a3 = antable["C0"] * SGsep_array ** antable["C1"] * API_array ** antable[
            "C2"] * Tres_array ** antable["C3"] * (
                     Pb_Velarde(API, SGsep, Tres, Psep, Tsep, Rsp) - 14.7) ** antable[
                 "C4"]

        Rsr = (a1 * Pr ** a2) + ((1 - a1) * Pr ** a3)
        Rs = Solution_GOR_Pb_ValkoMcCain(Psep, Tsep, API, Rsp) * Rsr

    elif P_array >= Pb_Velarde(API, SGsep, Tres, Psep, Tsep, Rsp):
        Rs = 0

    else:
        Rs = "Missing Data"

    if Units == 1:

        Rs = Rs

    else:
        Rs = Rs / (1.78107607 * 10 ** -1)

    return Rs


# print ("Rs =",Solution_GOR_Velarde2 (SGsep, API, Tres, P, Psep, Tsep, Rsp))

def Co_above_Pb_Spivey1(API, SGsep, P, Pres, Tres, Psep, Tsep, Rsp, Units=1):
    """The Spivey1 correlation calculates the Oil Compressibility from Pb to a higher
    pressure of interest, used for estimating certain fluid properties,the required
    input data is: Stock tank oil gravity [API], Separator gas Specific Gravity,
    Pressure [psia], Reservoir Pressure [psia] and Reservoir Temperature [F]. Also
    requires the Pb from the Velarde correlation and Rsb from the ValkoMcCain
    correlation

        If Units=1 use:
            Pressure in psia
            Temperature in F
            Stock-Tank oil gravity in API
            Solution Gas Oil Ratio @ Pb in Scf/STB

        If Units=2 use:
            Pressure in MPa
            Temperature in C
            Stock-Tank oil gravity in API
            Solution Gas Oil Ratio @ Pb in m3/m3"""

    if Units == 1:

        API_array = np.array(API)
        SGsep_array = np.array(SGsep)
        P_array = np.array(P)
        Tres_array = np.array(Tres)

    else:
        Tres_array = np.array(Tres)
        Tres_array = (Tres_array * 9 / 5) + 32
        API_array = np.array(API)
        SGsep_array = np.array(SGsep)
        P_array = np.array(P)
        P_array = P_array / (1.45037738 * 10 ** 2)

    if API_array > 0 and SGsep_array > 0 and Tres_array > 0 and P_array > 0:

        Pr = (P_array - 14.7) / (Pb_Velarde(API, SGsep, Tres, Psep, Tsep, Rsp) - 14.7)

        lnAPI = np.log(API_array)
        lnSGsep = np.log(SGsep_array)
        lnPb = np.log(Pb_Velarde(API, SGsep, Tres, Psep, Tsep, Rsp))
        lnPr = np.log(Pr)
        lnRsb = np.log(Solution_GOR_Pb_ValkoMcCain(Psep, Tsep, API, Rsp))
        lnTres = np.log(Tres_array)

        #         VARn            C0n                C1n                C2n
        Zntable = {"VAR1": lnAPI, "C01": 3.011, "C11": -2.6254, "C21": 0.497,
                   "VAR2": lnSGsep, "C02": -0.0835, "C12": -0.259, "C22": 0.382,
                   "VAR3": lnPb, "C03": 3.51, "C13": -0.0289, "C23": -0.0584,
                   "VAR4": lnPr, "C04": 0.327, "C14": -0.608, "C24": 0.0911,
                   "VAR5": lnRsb, "C05": -1.918, "C15": -0.642, "C25": 0.154,
                   "VAR6": lnTres, "C06": 2.52, "C16": -2.73, "C26": 0.429}

        z1 = Zntable["C01"] + (Zntable["C11"] * Zntable["VAR1"]) + (
                    Zntable["C21"] * Zntable["VAR1"] ** 2)
        z2 = Zntable["C02"] + (Zntable["C12"] * Zntable["VAR2"]) + (
                    Zntable["C22"] * Zntable["VAR2"] ** 2)
        z3 = Zntable["C03"] + (Zntable["C13"] * Zntable["VAR3"]) + (
                    Zntable["C23"] * Zntable["VAR3"] ** 2)
        z4 = Zntable["C04"] + (Zntable["C14"] * Zntable["VAR4"]) + (
                    Zntable["C24"] * Zntable["VAR4"] ** 2)
        z5 = Zntable["C05"] + (Zntable["C15"] * Zntable["VAR5"]) + (
                    Zntable["C25"] * Zntable["VAR5"] ** 2)
        z6 = Zntable["C06"] + (Zntable["C16"] * Zntable["VAR6"]) + (
                    Zntable["C26"] * Zntable["VAR6"] ** 2)

        sumZ = z1 + z2 + z3 + z4 + z5 + z6

        Cofb = np.exp(2.434 + (0.475 * sumZ) + (0.048 * sumZ ** 2) - np.log(10 ** 6))

    else:
        Cofb = "Missing Data"

    if Units == 1:
        Cofb = Cofb

    else:
        Cofb = Cofb / (1.45037738 * 10 ** 2)

    return Cofb


# print ("Cofb =",Co_above_Pb_Spivey1 (API, SGsep, P, Pres, Tres, Psep, Tsep, Rsp))

def Co_above_Pb_Spivey2(Pres, P, API, SGsep, Tres, Psep, Tsep, Rsp, Units=1):
    """The Spivey2 correlation calculates the Oil Compressibility from the initial
    reservoir pressure to a lower pressure but higher than Pb, used for material
    balance, the required imput data is: Pressure [psia], Reservoir Pressure [psia].
    Also requires the Pb from the Velarde correlation and Cofb and Cofbi from the
    Spivey1 correlation

        If Units=1 use:
            Pressure in psia
            Temperature in F
            Stock-Tank oil gravity in API
            Solution Gas Oil Ratio @ Pb in Scf/STB

        If Units=2 use:
            Pressure in MPa
            Temperature in C
            Stock-Tank oil gravity in API
            Solution Gas Oil Ratio @ Pb in m3/m3"""

    if Units == 1:

        Pres_array = np.array(Pres)
        P_array = np.array(P)
        API_array = np.array(API)
        SGsep_array = np.array(SGsep)
        Tres_array = np.array(Tres)

    else:
        Tres_array = np.array(Tres)
        Tres_array = (Tres_array * 9 / 5) + 32
        API_array = np.array(API)
        SGsep_array = np.array(SGsep)
        P_array = np.array(P)
        P_array = P_array / (1.45037738 * 10 ** 2)
        Pres_array = np.array(Pres)
        Pres_array = Pres_array / (1.45037738 * 10 ** 2)

    if Pres_array > 0 and P_array > 0:

        Pri = (Pres_array - 14.7) / (
                    Pb_Velarde(API, SGsep, Tres, Psep, Tsep, Rsp) - 14.7)

        lnAPI = np.log(API_array)
        lnSGsep = np.log(SGsep_array)
        lnPb = np.log(Pb_Velarde(API, SGsep, Tres, Psep, Tsep, Rsp))
        lnPri = np.log(Pri)
        lnRsb = np.log(Solution_GOR_Pb_ValkoMcCain(Psep, Tsep, API, Rsp))
        lnTres = np.log(Tres_array)
        #      VARn             C0n               C1n            C2n
        Zntable = {"VAR1": lnAPI, "C01": 3.011, "C11": -2.6254, "C21": 0.497,
                   "VAR2": lnSGsep, "C02": -0.0835, "C12": -0.259, "C22": 0.382,
                   "VAR3": lnPb, "C03": 3.51, "C13": -0.0289, "C23": -0.0584,
                   "VAR4": lnPri, "C04": 0.327, "C14": -0.608, "C24": 0.0911,
                   "VAR5": lnRsb, "C05": -1.918, "C15": -0.642, "C25": 0.154,
                   "VAR6": lnTres, "C06": 2.52, "C16": -2.73, "C26": 0.429}

        z1 = Zntable["C01"] + (Zntable["C11"] * Zntable["VAR1"]) + (
                    Zntable["C21"] * Zntable["VAR1"] ** 2)
        z2 = Zntable["C02"] + (Zntable["C12"] * Zntable["VAR2"]) + (
                    Zntable["C22"] * Zntable["VAR2"] ** 2)
        z3 = Zntable["C03"] + (Zntable["C13"] * Zntable["VAR3"]) + (
                    Zntable["C23"] * Zntable["VAR3"] ** 2)
        z4 = Zntable["C04"] + (Zntable["C14"] * Zntable["VAR4"]) + (
                    Zntable["C24"] * Zntable["VAR4"] ** 2)
        z5 = Zntable["C05"] + (Zntable["C15"] * Zntable["VAR5"]) + (
                    Zntable["C25"] * Zntable["VAR5"] ** 2)
        z6 = Zntable["C06"] + (Zntable["C16"] * Zntable["VAR6"]) + (
                    Zntable["C26"] * Zntable["VAR6"] ** 2)

        sumZ = z1 + z2 + z3 + z4 + z5 + z6

        Cofbi = np.exp(2.434 + (0.475 * sumZ) + (0.048 * sumZ ** 2) - np.log(10 ** 6))

        Cofi = (((Pb_Velarde(API, SGsep, Tres, Psep, Tsep,
                             Rsp) - P_array) * Co_above_Pb_Spivey1(API, SGsep, P, Pres,
                                                                   Tres, Psep, Tsep,
                                                                   Rsp)) - (
                        (Pb_Velarde(API, SGsep, Tres, Psep, Tsep,
                                    Rsp) - Pres_array) * Cofbi)) / (
                           Pres_array - P_array)

    else:
        Cofi = "Missing Data"

    if Units == 1:
        Cofi = Cofi

    else:
        Cofi = Cofi / (1.45037728 * 10 ** 2)

    return Cofi


# print ("Cofi =",Co_above_Pb_Spivey2 (Pres, P, API, SGsep, Tres, Psep, Tsep, Rsp))

def Co_above_Pb_Spivey3(P, API, SGsep, Tres, Psep, Tsep, Rsp, Pres, Units=1):
    """The Spivey3 correlation calculates the Oil Compressibility determined with the
    slope of a tangent line to the isothermal density-pressure curve at the pressure of
    interest, used for single phase pressure transient analysis, the required imput data
    is:
    Pressure [psia]. Also requires the Pb from the Velarde correlation, Cofb and Z
    constant values from the Spivey1 correlation

        If Units=1 use:
            Pressure in psia
            Temperature in F
            Stock-Tank oil gravity in API
            Solution Gas Oil Ratio @ Pb in Scf/STB

        If Units=2 use:
            Pressure in MPa
            Temperature in C
            Stock-Tank oil gravity in API
            Solution Gas Oil Ratio @ Pb in m3/m3"""

    if Units == 1:

        P_array = np.array(P)
        API_array = np.array(API)
        SGsep_array = np.array(SGsep)
        Tres_array = np.array(Tres)

    else:
        Tres_array = np.array(Tres)
        Tres_array = (Tres_array * 9 / 5) + 32
        API_array = np.array(API)
        SGsep_array = np.array(SGsep)
        P_array = np.array(P)
        P_array = P_array / (1.45037738 * 10 ** 2)

    if P_array > 0:

        Pr = (P_array - 14.7) / (Pb_Velarde(API, SGsep, Tres, Psep, Tsep, Rsp) - 14.7)

        dZdP = (-0.608 + (0.1822 * np.log(
            P_array / Pb_Velarde(API, SGsep, Tres, Psep, Tsep, Rsp)))) / P_array

        lnAPI = np.log(API_array)
        lnSGsep = np.log(SGsep_array)
        lnPb = np.log(Pb_Velarde(API, SGsep, Tres, Psep, Tsep, Rsp))
        lnPr = np.log(Pr)
        lnRsb = np.log(Solution_GOR_Pb_ValkoMcCain(Psep, Tsep, API, Rsp))
        lnTres = np.log(Tres_array)

        #         VARn            C0n                C1n                C2n
        Zntable = {"VAR1": lnAPI, "C01": 3.011, "C11": -2.6254, "C21": 0.497,
                   "VAR2": lnSGsep, "C02": -0.0835, "C12": -0.259, "C22": 0.382,
                   "VAR3": lnPb, "C03": 3.51, "C13": -0.0289, "C23": -0.0584,
                   "VAR4": lnPr, "C04": 0.327, "C14": -0.608, "C24": 0.0911,
                   "VAR5": lnRsb, "C05": -1.918, "C15": -0.642, "C25": 0.154,
                   "VAR6": lnTres, "C06": 2.52, "C16": -2.73, "C26": 0.429}

        z1 = Zntable["C01"] + (Zntable["C11"] * Zntable["VAR1"]) + (
                    Zntable["C21"] * Zntable["VAR1"] ** 2)
        z2 = Zntable["C02"] + (Zntable["C12"] * Zntable["VAR2"]) + (
                    Zntable["C22"] * Zntable["VAR2"] ** 2)
        z3 = Zntable["C03"] + (Zntable["C13"] * Zntable["VAR3"]) + (
                    Zntable["C23"] * Zntable["VAR3"] ** 2)
        z4 = Zntable["C04"] + (Zntable["C14"] * Zntable["VAR4"]) + (
                    Zntable["C24"] * Zntable["VAR4"] ** 2)
        z5 = Zntable["C05"] + (Zntable["C15"] * Zntable["VAR5"]) + (
                    Zntable["C25"] * Zntable["VAR5"] ** 2)
        z6 = Zntable["C06"] + (Zntable["C16"] * Zntable["VAR6"]) + (
                    Zntable["C26"] * Zntable["VAR6"] ** 2)

        sumZ = z1 + z2 + z3 + z4 + z5 + z6

        dCofbdP = Co_above_Pb_Spivey1(API, SGsep, P, Pres, Tres, Psep, Tsep, Rsp) * (
                    (0.475 + (0.096 * sumZ)) * dZdP)
        Co = Co_above_Pb_Spivey1(API, SGsep, P, Pres, Tres, Psep, Tsep, Rsp) + (
                    (P_array - Pb_Velarde(API, SGsep, Tres, Psep, Tsep, Rsp)) * dCofbdP)

    else:
        Co = "Missing Data"

    if Units == 1:
        Co = Co

    else:
        Co = Co / (1.45037728 * 10 ** 2)

    return Co


# print ("Co =",Co_above_Pb_Spivey3 (P, API, SGsep, Tres, Psep, Tsep, Rsp, Pres))

# ----CONTINUE AFTER FVF CORRELATIONS----

def Den_Pb(SGsep, Tres, Psep, Tsep, API, Rsp, Units=1):
    """The Den_Pb function calculates the reservoir oil Densities at reservoir
    pressures equal to bubble point pressure, the required imput data is: Separator gas
    Specific Gravity and Reservoir Temperature [F]. Also requires the Rsb from
    ValkoMcCain correlation, SG and SGst from ValkoMcCain2 correlation and Pb from the
    Velarde correlation

        If Units=1 use:
            Pressure in psia
            Temperature in F
            Stock-Tank oil gravity in API
            Solution Gas Oil Ratio @ Pb in Scf/STB

        If Units=2 use:
            Pressure in MPa
            Temperature in C
            Stock-Tank oil gravity in API
            Solution Gas Oil Ratio @ Pb in m3/m3"""

    if Units == 1:

        SGsep_array = np.array(SGsep)
        Tres_array = np.array(Tres)

    else:
        Tres_array = np.array(Tres)
        Tres_array = (Tres_array * 9 / 5) + 32
        SGsep_array = np.array(SGsep)

    Den_pseudo1 = 52.8 - (0.01 * (Solution_GOR_Pb_ValkoMcCain(Psep, Tsep, API, Rsp)))

    if SGsep_array > 0 and Tres_array > 0:

        antable = {"a0": -49.8930, "a1": 85.0149, "a2": -3.70373, "a3": 0.0479818,
                   "a4": 2.98914, "a5": -0.0356888}

        den_ap = antable["a0"] + antable["a1"] * SGsep_array + antable[
            "a2"] * SGsep_array * Den_pseudo1 + antable[
                     "a3"] * SGsep_array * Den_pseudo1 ** 2 + antable[
                     "a4"] * Den_pseudo1 + antable["a5"] * Den_pseudo1 ** 2

        Den_pseudo2 = ((Solution_GOR_Pb_ValkoMcCain(Psep, Tsep, API,
                                                    Rsp) * Spec_grav_ValkoMcCain2(Psep,
                                                                                  Tsep,
                                                                                  API,
                                                                                  Rsp,
                                                                                  SGsep)) + (
                               4600 * Spec_grav_st_ValkoMcCain2(Psep, Tsep, API, Rsp,
                                                                SGsep))) / (73.71 + (
                (Solution_GOR_Pb_ValkoMcCain(Psep, Tsep, API,
                                             Rsp) * Spec_grav_ValkoMcCain2(Psep, Tsep,
                                                                           API, Rsp,
                                                                           SGsep)) / den_ap))

        while ((abs(Den_pseudo1 - Den_pseudo2)) / Den_pseudo2) > 0.1:

            Den_pseudo1 = Den_pseudo2

            den_ap = antable["a0"] + antable["a1"] * SGsep_array + antable[
                "a2"] * SGsep_array * Den_pseudo1 + antable[
                         "a3"] * SGsep_array * Den_pseudo1 ** 2 + antable[
                         "a4"] * Den_pseudo1 + antable["a5"] * Den_pseudo1 ** 2

            Den_pseudo2 = ((Solution_GOR_Pb_ValkoMcCain(Psep, Tsep, API,
                                                        Rsp) * Spec_grav_ValkoMcCain2(
                Psep, Tsep, API, Rsp, SGsep)) + (
                                   4600 * Spec_grav_st_ValkoMcCain2(Psep, Tsep, API,
                                                                    Rsp, SGsep))) / (
                                      73.71 + (
                                      (Solution_GOR_Pb_ValkoMcCain(Psep, Tsep, API,
                                                                   Rsp) * Spec_grav_ValkoMcCain2(
                                          Psep, Tsep, API, Rsp, SGsep)) / den_ap))

            if ((abs(Den_pseudo1 - Den_pseudo2)) / Den_pseudo2) < 0.1:
                break

        den_adjP = ((0.167 + 16.181 * (10 ** (-0.0425 * Den_pseudo2))) * (
                    Pb_Velarde(API, SGsep, Tres, Psep, Tsep, Rsp) / 1000)) - 0.01 * (
                           (0.299 + 263 * (10 ** (-0.0603 * Den_pseudo2))) * (
                               Pb_Velarde(API, SGsep, Tres, Psep, Tsep,
                                          Rsp) / 1000) ** 2)

        den_fake = den_adjP + Den_pseudo2

        den_adjT = ((0.00302 + (1.505 * den_fake ** -0.951)) * (
                    (Tres_array - 60) ** 0.938)) - (
                           (0.0216 - 0.0233 * (10 ** (-0.0161 * den_fake))) * (
                               (Tres_array - 60) ** 0.475))

        den_Pb = den_fake - den_adjT

    if Units == 1:
        den_Pb = den_Pb

    else:
        den_Pb = den_Pb / (1.60184634 * 10 ** -2)

    return den_Pb


# print ("den_Pb =",Den_Pb (SGsep, Tres, Psep, Tsep, API, Rsp))

def Den_underPb(SGsep, Tres, P, API, Psep, Tsep, Rsp, Units=1):
    """The Den_underPb function calculates the reservoir oil Densities at reservoir
    pressures less than bubble point pressure, the required imput data is: Separator
    gas Specific Gravity, Reservoir Temperature [F] and Pressure [psia]. Also requires
    the Rs from Velarde2 correlation and SG and SGst from ValkoMcCain2 correlation

        If Units=1 use:
            Pressure in psia
            Temperature in F
            Stock-Tank oil gravity in API
            Solution Gas Oil Ratio @ Pb in Scf/STB

        If Units=2 use:
            Pressure in MPa
            Temperature in C
            Stock-Tank oil gravity in API
            Solution Gas Oil Ratio @ Pb in m3/m3"""

    if Units == 1:
        SGsep_array = np.array(SGsep)
        Tres_array = np.array(Tres)
        P_array = np.array(P)

    else:
        Tres_array = np.array(Tres)
        Tres_array = (Tres_array * 9 / 5) + 32
        SGsep_array = np.array(SGsep)
        P_array = np.array(P)
        P_array = P_array / (1.45037738 * 10 ** 2)

    Den_pseudo1 = 52.8 - (
                0.01 * (Solution_GOR_Velarde2(SGsep, API, Tres, P, Psep, Tsep, Rsp)))

    if SGsep_array > 0 and Tres_array > 0:

        antable = {"a0": -49.8930, "a1": 85.0149, "a2": -3.70373, "a3": 0.0479818,
                   "a4": 2.98914, "a5": -0.0356888}

        den_ap = antable["a0"] + antable["a1"] * SGsep_array + antable[
            "a2"] * SGsep_array * Den_pseudo1 + antable[
                     "a3"] * SGsep_array * Den_pseudo1 ** 2 + antable[
                     "a4"] * Den_pseudo1 + antable["a5"] * Den_pseudo1 ** 2

        Den_pseudo2 = ((Solution_GOR_Velarde2(SGsep, API, Tres, P, Psep, Tsep,
                                              Rsp) * Spec_grav_ValkoMcCain2(Psep, Tsep,
                                                                            API, Rsp,
                                                                            SGsep)) + (
                                   4600 * Spec_grav_st_ValkoMcCain2(Psep, Tsep, API,
                                                                    Rsp, SGsep))) / (
                              73.71 + ((Solution_GOR_Velarde2(SGsep, API, Tres, P, Psep,
                                                              Tsep,
                                                              Rsp) * Spec_grav_ValkoMcCain2(
                          Psep, Tsep, API, Rsp, SGsep)) / den_ap))

        while ((abs(Den_pseudo1 - Den_pseudo2)) / Den_pseudo2) * 100 > 10:

            Den_pseudo1 = Den_pseudo2

            den_ap = antable["a0"] + antable["a1"] * SGsep_array + antable[
                "a2"] * SGsep_array * Den_pseudo1 + antable[
                         "a3"] * SGsep_array * Den_pseudo1 ** 2 + antable[
                         "a4"] * Den_pseudo1 + antable["a5"] * Den_pseudo1 ** 2

            Den_pseudo2 = ((Solution_GOR_Velarde2(SGsep, API, Tres, P, Psep, Tsep,
                                                  Rsp) * Spec_grav_ValkoMcCain2(Psep,
                                                                                Tsep,
                                                                                API,
                                                                                Rsp,
                                                                                SGsep)) + (
                                       4600 * Spec_grav_st_ValkoMcCain2(Psep, Tsep, API,
                                                                        Rsp,
                                                                        SGsep))) / (
                                  73.71 + ((Solution_GOR_Velarde2(SGsep, API, Tres, P,
                                                                  Psep, Tsep,
                                                                  Rsp) * Spec_grav_ValkoMcCain2(
                              Psep, Tsep, API, Rsp, SGsep)) / den_ap))

            if ((abs(Den_pseudo1 - Den_pseudo2)) / Den_pseudo2) < 0.1:
                break

        den_adjP = ((0.167 + 16.181 * (10 ** (-0.0425 * Den_pseudo2))) * (
                    P_array / 1000)) - 0.01 * (
                           (0.299 + 263 * (10 ** (-0.0603 * Den_pseudo2))) * (
                               P_array / 1000) ** 2)

        den_fake = den_adjP + Den_pseudo2

        den_adjT = ((0.00302 + (1.505 * den_fake ** -0.951)) * (
                    (Tres_array - 60) ** 0.938)) - (
                           (0.0216 - 0.0233 * (10 ** (-0.0161 * den_fake))) * (
                               (Tres_array - 60) ** 0.475))

        den_underPb = den_fake - den_adjT

    if Units == 1:
        den_underPb = den_underPb

    else:
        den_underPb = den_underPb / (1.60184634 * 10 ** -2)

    return den_underPb


# print ("den_underPb =",Den_underPb (SGsep, Tres, P, API, Psep, Tsep, Rsp))

def Den_abovePb(P, SGsep, Tres, Psep, Tsep, API, Rsp, Pres, Units=1):
    """The Den_abovePb function calculates the reservoir oil Densities at reservoir
    pressures greater than bubble point pressure, the required imput data is: Pressure
    [psia]. Also requires the Den_Pb function, Cofb from Spivey1 correlation and Pb from
    Velarde correlation

        If Units=1 use:
            Pressure in psia
            Temperature in F
            Stock-Tank oil gravity in API
            Solution Gas Oil Ratio @ Pb in Scf/STB

        If Units=2 use:
            Pressure in MPa
            Temperature in C
            Stock-Tank oil gravity in API
            Solution Gas Oil Ratio @ Pb in m3/m3"""

    if Units == 1:

        P_array = np.array(P)

    else:

        P_array = np.array(P)
        P_array = P_array / (1.45037738 * 10 ** 2)

    den_abovePb = Den_Pb(SGsep, Tres, Psep, Tsep, API, Rsp) * np.exp(
        Co_above_Pb_Spivey1(API, SGsep, P, Pres, Tres, Psep, Tsep, Rsp) * (
                    P_array - Pb_Velarde(API, SGsep, Tres, Psep, Tsep, Rsp)))

    if Units == 1:
        den_abovePb = den_abovePb

    else:
        den_abovePb = den_abovePb / (1.60184634 * 10 ** -2)

    return den_abovePb


# print("den_abovePb =",Den_abovePb(P, SGsep, Tres, Psep, Tsep, API, Rsp, Pres))

def FVF_Pb(Den_sto, Psep, Tsep, API, Rsp, SGsep, Tres, Units=1):
    """The FVF_Pb function calculates the Formation Volume Factor at reservoir pressure
    equal to bubblepoint pressure, the required imput data is: Density of stock-tank
    oil [lb/cu ft]. Also requires the Rsb from ValkoMcCain correlation, SG from
    ValkoMcCain2 correlation, Den_Pb function

        If Units=1 use:
            Pressure in psia
            Temperature in F
            Stock-Tank oil gravity in API
            Solution Gas Oil Ratio @ Pb in Scf/STB
            Density @ Stock-Tank in lb/ft3

        If Units=2 use:
            Pressure in MPa
            Temperature in C
            Stock-Tank oil gravity in API
            Solution Gas Oil Ratio @ Pb in m3/m3
            Density @ Stock-Tank in gr/cm3"""

    if Units == 1:

        Den_sto_array = np.array(Den_sto)

    else:
        Den_sto_array = np.array(Den_sto)
        Den_sto_array = Den_sto_array / (1.60184634 * 10 ** -2)

    Bob = (Den_sto_array + (
            0.01357 * Solution_GOR_Pb_ValkoMcCain(Psep, Tsep, API,
                                                  Rsp) * Spec_grav_st_ValkoMcCain2(Psep,
                                                                                   Tsep,
                                                                                   API,
                                                                                   Rsp,
                                                                                   SGsep))) / Den_Pb(
        SGsep, Tres, Psep, Tsep, API, Rsp)

    return Bob


# print ("Bob =",FVF_Pb (Den_sto, Psep, Tsep, API, Rsp, SGsep, Tres))

def FVF_underPb(Den_sto, Psep, Tsep, API, Rsp, SGsep, Tres, P, Units=1):
    """The FVF_underPb function calculates the Formation Volume Factor at reservoir
    pressure less than bubblepoint pressure, the required imput data is: Density of
    stock-tank oil [lb/cu ft]. Also requires the Rs from Velarde correlation, SG from
    ValkoMcCain2 correlation and Den_underPb function

        If Units=1 use:
            Pressure in psia
            Temperature in F
            Stock-Tank oil gravity in API
            Solution Gas Oil Ratio @ Pb in Scf/STB
            Density @ Stock-Tank in lb/ft3

        If Units=2 use:
            Pressure in MPa
            Temperature in C
            Stock-Tank oil gravity in API
            Solution Gas Oil Ratio @ Pb in m3/m3
            Density @ Stock-Tank in gr/cm3"""

    if Units == 1:

        Den_sto_array = np.array(Den_sto)

    else:
        Den_sto_array = np.array(Den_sto)
        Den_sto_array = Den_sto_array / (1.60184634 * 10 ** -2)

    Bo_underPb = (Den_sto_array + (
            0.01357 * Solution_GOR_Velarde2(SGsep, API, Tres, P, Psep, Tsep,
                                            Rsp) * Spec_grav_st_ValkoMcCain2(Psep, Tsep,
                                                                             API, Rsp,
                                                                             SGsep))) / Den_underPb(
        SGsep, Tres, P, API, Psep, Tsep, Rsp)

    return Bo_underPb


# print ("Bo_underPb =",FVF_underPb(Den_sto, Psep, Tsep, API, Rsp, SGsep, Tres, P))

def FVF_abovePb(P, Den_sto, Psep, Tsep, API, Rsp, SGsep, Tres, Pres, Units=1):
    """The FVF_underPb function calculates the Formation Volume Factor at reservoir
    pressure greater than bubblepoint pressure, the required imput data is: Pressure
    [psia]. Also requires the Bob from FVF_Pb function, Cofb from Spivey1 correlation
    and Pb from Velarde correlation

        If Units=1 use:
            Pressure in psia
            Temperature in F
            Stock-Tank oil gravity in API
            Solution Gas Oil Ratio @ Pb in Scf/STB
            Density @ Stock-Tank in lb/ft3

        If Units=2 use:
            Pressure in MPa
            Temperature in C
            Stock-Tank oil gravity in API
            Solution Gas Oil Ratio @ Pb in m3/m3
            Density @ Stock-Tank in gr/cm3"""

    if Units == 1:

        P_array = np.array(P)

    else:
        P_array = np.array(P)
        P_array = P_array / (1.45037738 * 10 ** 2)

    Bo_abovePb = FVF_Pb(Den_sto, Psep, Tsep, API, Rsp, SGsep, Tres) * np.exp(
        Co_above_Pb_Spivey1(API, SGsep, P, Pres, Tres, Psep, Tsep, Rsp) * (
                    Pb_Velarde(API, SGsep, Tres, Psep, Tsep, Rsp) - P_array))

    return Bo_abovePb


# print ("Bo_abovePb =",FVF_abovePb (P, Den_sto, Psep, Tsep, API, Rsp, SGsep, Tres, Pres))

# ----CONTINUE OF VISCOSITY CORRELATIONS----


def Viscosity_Pb(API, Tres, Psep, Tsep, Rsp, Units=1):
    """The Viscosity_Pb function calculates the Viscosity at reservoir pressures equal
     to bubblepoint pressure, the required imput data is: Stock tank oil gravity [API]
      and Reservoir Temperature [F]. Also requires the Rsb from ValkoMcCain correlation

        If Units=1 use:
            Pressure in psia
            Temperature in F
            Stock-Tank oil gravity in API
            Solution Gas Oil Ratio @ Pb in Scf/STB

        If Units=2 use:
            Pressure in MPa
            Temperature in C
            Stock-Tank oil gravity in API
            Solution Gas Oil Ratio @ Pb in m3/m3"""

    if Units == 1:

        API_array = np.array(API)
        Tres_array = np.array(Tres)

    else:
        API_array = np.array(API)
        Tres_array = np.array(Tres)
        Tres_array = (Tres_array * 9 / 5) + 32

    C = (10 ** (3.0324 - 0.02023 * API_array)) * (Tres_array ** -1.163)

    uoD = (10 ** C) - 1

    A = 10.715 * (Solution_GOR_Pb_ValkoMcCain(Psep, Tsep, API, Rsp) + 100) ** -0.515
    B = 5.44 * (Solution_GOR_Pb_ValkoMcCain(Psep, Tsep, API, Rsp) + 150) ** -0.338

    uob = A * uoD ** B

    if Units == 1:
        uob = uob

    else:
        uob = uob / (0.001)

    return uob


# print("uob =",Viscosity_Pb (API, Tres, Psep, Tsep, Rsp))

def Viscosity_underPb(SGsep, API, Tres, P, Psep, Tsep, Rsp, Units=1):
    """The Viscosity_underPb function calculates the Viscosity at reservoir pressures
    under bubblepoint pressure. Requires the Rs from Velarde2 correlation and uoR from
     Viscosity_Pb function

        If Units=1 use:
            Pressure in psia
            Temperature in F
            Stock-Tank oil gravity in API
            Solution Gas Oil Ratio @ Pb in Scf/STB

        If Units=2 use:
            Pressure in MPa
            Temperature in C
            Stock-Tank oil gravity in API
            Solution Gas Oil Ratio @ Pb in m3/m3"""

    if Units == 1:

        API_array = np.array(API)
        Tres_array = np.array(Tres)

    else:
        API_array = np.array(API)
        Tres_array = np.array(Tres)
        Tres_array = (Tres_array * 9 / 5) + 32

    A = 10.715 * (Solution_GOR_Velarde2(SGsep, API, Tres, P, Psep, Tsep,
                                        Rsp) + 100) ** -0.515
    B = 5.44 * (Solution_GOR_Velarde2(SGsep, API, Tres, P, Psep, Tsep,
                                      Rsp) + 150) ** -0.338

    C = (10 ** (3.0324 - 0.02023 * API_array)) * (Tres_array ** -1.163)

    uoD = (10 ** C) - 1

    uoR = A * uoD ** B

    if Units == 1:
        uoR = uoR

    else:
        uoR = uoR / (0.001)

    return uoR


# print("uoR_underPb =",Viscosity_underPb (SGsep, API, Tres, P, Psep, Tsep, Rsp))

def Viscosity_abovePb(P, API, Tres, Psep, Tsep, Rsp, SGsep, Units=1):
    """The Viscosity_abovePb function calculates the Viscosity at reservoir pressures
    above bubblepoint pressure, the required imput data is: Pressure [psia]. Also
    requires the uob from Viscosity_Pb function and Pb from Velarde correlation

        If Units=1 use:
            Pressure in psia
            Temperature in F
            Stock-Tank oil gravity in API
            Solution Gas Oil Ratio @ Pb in Scf/STB

        If Units=2 use:
            Pressure in MPa
            Temperature in C
            Stock-Tank oil gravity in API
            Solution Gas Oil Ratio @ Pb in m3/m3"""

    if Units == 1:

        P_array = np.array(P)

    else:
        P_array = np.array(P)
        P_array = P_array / (1.45037738 * 10 ** 2)

    A = -1.0146 + (1.3322 * np.log10(Viscosity_Pb(API, Tres, Psep, Tsep, Rsp))) - (
                0.4876 * (np.log10(Viscosity_Pb(API, Tres, Psep, Tsep, Rsp))) ** 2) - (
                1.15036 * (np.log10(Viscosity_Pb(API, Tres, Psep, Tsep, Rsp))) ** 3)

    uoR = Viscosity_Pb(API, Tres, Psep, Tsep, Rsp) + ((1.3449 * 10 ** -3) * (
                P_array - Pb_Velarde(API, SGsep, Tres, Psep, Tsep, Rsp)) * 10 ** A)

    if Units == 1:
        uoR = uoR

    else:
        uoR = uoR / (0.001)

    return uoR
# print("uoR_abovePb =",Viscosity_abovePb (P, API, Tres, Psep, Tsep, Rsp, SGsep))


"""----------------------------Water Correlations------------------------------------"""

def RS_bw(P, T, salinity, unit=1):
    """This function calculate the Rs of methane dissolved brine
    The user must select UNITS 1 to English System and 2 to International System

    IF UNITS 1 is selected introduce:
        | Pressure (P) in PSI
        | Temperature (T) in ºF
        | Salinity in Weight Fraction
        | Rs returns in scf/STB

    IF UNITS 2 is selected introduce:
        | Pressure (P) in MPa
        | Temperature (T) in ºC
        | Salinity in Weight Fraction
        | Rs returns in cm3/cm3"""

    if unit == 1:
        P_array = np.array(P)
        P_array = P_array * (6.89475729 * 10 ** -3.0)
        T_array = np.array(T)
        T_array = (T_array - 32) * 5 / 9
        m = np.array(salinity)
        m = 1000 * (m / 1000000) / (58.44 * (1 - (m / 1000000)))

    else:
        P_array = np.array(P)
        T_array = np.array(T)
        m = np.array(salinity)
        m = 1000 * (m / 1000000) / (58.44 * (1 - (m / 1000000)))

    # Step1: Calculate the vapor pressure of pure water

    T_array = T_array + 273

    # Coefficients Table 4_9
    a1Po = -7.85951783
    a2Po = 1.84408259
    a3Po = -11.7866497
    a4Po = 22.6807411
    a5Po = -15.9618719
    a6Po = 1.80122502

    # Ec 4.13
    vad = 1 - T_array / 674.096  # Adimensional T EC.4.14
    # Ec 4.14
    Po = 22.064 * (np.exp((674.096 / T_array) * (
                a1Po * vad + a2Po * vad ** 1.5 + a3Po * vad ** 3 + a4Po * vad ** 3.5 + a5Po * vad ** 4 + a6Po * vad ** 7.5)))

    # Step2: Calculate the solubility of methane in pure water at (P,T)

    T_array = T_array - 273

    # Coefficients Table 4_10 EC 4.1
    a1A_T = 0
    a2A_T = -0.004462
    a3A_T = -0.06763
    a4A_T = 0
    a5A_T = 0

    A_T = (a1A_T * (T_array / 100) ** 2 + a2A_T * (T_array / 100) + a3A_T) / (
                a4A_T * (T_array / 100) ** 2 + a5A_T * (T_array / 100) + 1)

    a1B_T = -0.03602
    a2B_T = 0.18917
    a3B_T = 0.97242
    a4B_T = 0
    a5B_T = 0

    B_T = (a1B_T * (T_array / 100) ** 2 + a2B_T * (T_array / 100) + a3B_T) / (
                a4B_T * (T_array / 100) ** 2 + a5B_T * (T_array / 100) + 1)

    a1C_T = 0.6855
    a2C_T = -3.1992
    a3C_T = -3.7968
    a4C_T = 0.07711
    a5C_T = 0.2229

    C_T = (a1C_T * (T_array / 100) ** 2 + a2C_T * (T_array / 100) + a3C_T) / (
                a4C_T * (T_array / 100) ** 2 + a5C_T * (T_array / 100) + 1)

    # Ec 4.13 (Solubility of methane in pure water)
    mCH4pw = np.exp(A_T * (np.log(P_array - Po)) ** 2 + B_T * np.log(P_array - Po) + C_T)  # Ec 4.15

    T_array = T_array + 273

    # Step3: Calculate the solubility of methane in brine water at (P,T)

    # Coefficients Table 4_11
    lambC1 = -0.80898
    lambC2 = 1.0827 * 10 ** (-3)
    lambC3 = 183.85
    lambC6 = 3.924 * 10 ** (-4)
    lambC10 = -1.97 * 10 ** (-6)

    # Ec 4.16
    lambCH4 = lambC1 + lambC2 * T_array + lambC3 / T_array + lambC6 * P_array + lambC10 * P_array ** 2  # EC4.16

    # Ec 4.18
    mCH4bw = mCH4pw * np.exp(-2 * lambCH4 * m - (-3.89 * 10 ** (-3)) * m ** 2)

    Z = 0.98 # TODO input for user?

    # Ec 4.34
    VMCH4gsc = Z * 8.314467 * (15.5 + 273) / 0.001

    # Step4: Calculate the formation volumen factor

    # Ec 4.36

    if unit == 1:
        Rs = (mCH4bw * VMCH4gsc) / ((1000 + m * 58.44) * 0.99681786)
        Rs = Rs * 5.61458333

    else:
        Rs = (mCH4bw * VMCH4gsc) / ((1000 + m * 58.44) * 0.99681786)

    return Rs


def Bo_bw(P, T, salinity, unit=1):
    """This function calculate the Bo of methane dissolved brine
    The user must select UNITS 1 to Field units and 2 to International System

    IF UNITS 1 is selected introduce:
        | Pressure (P) in PSI
        | Temperature (T) in ºF
        | Salinity in ppm
        | Bo returns in bbl/STB

    IF UNITS 2 is selected introduce:
        Pressure (P) in MPa
        Temperature (T) in ºC
        Salinity in ppm
        Bo returns in cm3/cm3"""

    # Step1: Calculate the derivatives of Duan

    if unit == 1:
        P_array = np.array(P)
        P_array = P_array * (6.89475729 * 10 ** -3.0)
        T_array = np.array(T)
        T_array = (T_array - 32) * 5 / 9
        m = np.array(salinity)
        m = 1000 * (m / 1000000) / (58.44 * (1 - (m / 1000000)))

    else:
        P_array = np.array(P)
        T_array = np.array(T)
        m = np.array(salinity)
        m = 1000 * (m / 1000000) / (58.44 * (1 - (m / 1000000)))

    T_array = T_array + 273
    # Coefficients Table 4_11
    UCH4c6 = 7.6985890 * 10 ** (-2)
    UCH4c7 = -5.0253331 * 10 ** (-5)
    UCH4c8 = -30.092013
    UCH4c9 = 4.8468502 * 10 ** (3)

    CH4c6 = 3.92 * 10 ** (-4)
    CH4c10 = -1.97 * 10 ** (-6)

    # Ec 4.19
    dch4dp = UCH4c6 + UCH4c7 * T_array + UCH4c8 / T_array + UCH4c9 / (T_array) ** 2
    # Ec 4.20
    dlamdp = CH4c6 + 2 * CH4c10 * P_array
    # Ec 4.22 Partial molar volume of methane in brine
    VMCH4B = 8.314467 * T_array * (dch4dp + 2 * m * dlamdp + m ** 2 * 0)

    # Step2: Calculate the density of methane-free brine water

    T_array = T_array - 273

    # Coefficients Table 4_6 and Ec 4.1
    a1_T_70 = -0.127213
    a2_T_70 = 0.645486
    a3_T_70 = 1.03265
    a4_T_70 = -0.070291
    a5_T_70 = 0.639589

    a_T_70 = (a1_T_70 * (T_array / 100) ** 2 + a2_T_70 * (T_array / 100) + a3_T_70) / (
                a4_T_70 * (T_array / 100) ** 2 + a5_T_70 * (T_array / 100) + 1)

    a1_EW = 4.221
    a2_EW = -3.478
    a3_EW = 6.221
    a4_EW = 0.5182
    a5_EW = -0.4405

    a_T_EW = (a1_EW * (T_array / 100) ** 2 + a2_EW * (T_array / 100) + a3_EW) / (
                a4_EW * (T_array / 100) ** 2 + a5_EW * (T_array / 100) + 1)

    a1_FW = -11.403
    a2_FW = 29.932
    a3_FW = 27.952
    a4_FW = 0.20684
    a5_FW = 0.3768

    a_T_FW = (a1_FW * (T_array / 100) ** 2 + a2_FW * (T_array / 100) + a3_FW) / (
                a4_FW * (T_array / 100) ** 2 + a5_FW * (T_array / 100) + 1)

    # Coefficients Table 4_7 and Ec 4.1

    a1_DM2 = -1.1149 * 10 ** (-4)
    a2_DM2 = 1.7105 * 10 ** (-4)
    a3_DM2 = -4.3766 * 10 ** (-4)
    a4_DM2 = 0
    a5_DM2 = 0

    a_T_DM2 = (a1_DM2 * (T_array / 100) ** 2 + a2_DM2 * (T_array / 100) + a3_DM2) / (
            a4_DM2 * (T_array / 100) ** 2 + a5_DM2 * (T_array / 100) + 1)

    a1_DM3_2 = -8.878 * 10 ** (-4)
    a2_DM3_2 = -1.388 * 10 ** (-4)
    a3_DM3_2 = -2.96318 * 10 ** (-3)
    a4_DM3_2 = 0
    a5_DM3_2 = 0.51103

    a_T_DM3_2 = (a1_DM3_2 * (T_array / 100) ** 2 + a2_DM3_2 * (T_array / 100) + a3_DM3_2) / (
                a4_DM3_2 * (T_array / 100) ** 2 + a5_DM3_2 * (T_array / 100) + 1)

    a1_DM1 = 2.1466 * 10 ** (-3)
    a2_DM1 = 1.2427 * 10 ** (-2)
    a3_DM1 = 4.2648 * 10 ** (-2)
    a4_DM1 = -8.1009 * 10 ** (-2)
    a5_DM1 = 0.525417

    a_T_DM1 = (a1_DM1 * (T_array / 100) ** 2 + a2_DM1 * (T_array / 100) + a3_DM1) / (
                a4_DM1 * (T_array / 100) ** 2 + a5_DM1 * (T_array / 100) + 1)

    a1_DM1_2 = 2.366 * 10 ** (-4)
    a2_DM1_2 = -3.636 * 10 ** (-4)
    a3_DM1_2 = -2.278 * 10 ** (-4)
    a4_DM1_2 = 0
    a5_DM1_2 = 0

    a_T_DM1_2 = (a1_DM1_2 * (T_array / 100) ** 2 + a2_DM1_2 * (T_array / 100) + a3_DM1_2) / (
                a4_DM1_2 * (T_array / 100) ** 2 + a5_DM1_2 * (T_array / 100) + 1)

    # Ec 4.6 Density of brine water at (70MPa and T)
    dbw_T_70 = a_T_70 + a_T_DM2 * 1 * m ** 2 + a_T_DM3_2 * 1 * m ** (3 / 2) + a_T_DM1 * 1 * m + a_T_DM1_2 * 1 * m ** (
                1 / 2)

    # Coefficients Table 4_8 Ec 4.1

    a1_EM = 0
    a2_EM = 0
    a3_EM = 0.1249
    a4_EM = 0
    a5_EM = 0

    a_T_EM = (a1_EM * (T_array / 100) ** 2 + a2_EM * (T_array / 100) + a3_EM) / (
            a4_EM * (T_array / 100) ** 2 + a5_EM * (T_array / 100) + 1)

    a1_FM3_2 = -0.617
    a2_FM3_2 = -0.747
    a3_FM3_2 = -0.4339
    a4_FM3_2 = 0
    a5_FM3_2 = 10.26

    a_T_FM3_2 = (a1_FM3_2 * (T_array / 100) ** 2 + a2_FM3_2 * (T_array / 100) + a3_FM3_2) / (
            a4_FM3_2 * (T_array / 100) ** 2 + a5_FM3_2 * (T_array / 100) + 1)

    a1_FM1 = 0
    a2_FM1 = 9.917
    a3_FM1 = 5.1128
    a4_FM1 = 0
    a5_FM1 = 3.892

    a_T_FM1 = (a1_FM1 * (T_array / 100) ** 2 + a2_FM1 * (T_array / 100) + a3_FM1) / (
                a4_FM1 * (T_array / 100) ** 2 + a5_FM1 * (T_array / 100) + 1)

    a1_FM1_2 = 0.0365
    a2_FM1_2 = -0.0369
    a3_FM1_2 = 0
    a4_FM1_2 = 0
    a5_FM1_2 = 0

    a_T_FM1_2 = (a1_FM1_2 * (T_array / 100) ** 2 + a2_FM1_2 * (T_array / 100) + a3_FM1_2) / (
                a4_FM1_2 * (T_array / 100) ** 2 + a5_FM1_2 * (T_array / 100) + 1)

    Eb_T_M = a_T_EW + a_T_EM * m

    FB_T_M = a_T_FW + a_T_FM3_2 * m ** (3 / 2) + a_T_FM1 * m + a_T_FM1_2 * m ** (1 / 2)

    Ib_T_70_m = (1 / Eb_T_M) * np.log(np.abs(Eb_T_M + FB_T_M))

    Ib_T_P_m = (1 / Eb_T_M) * np.log(np.abs(Eb_T_M * (P_array / 70) + FB_T_M))

    # Density of methane-free brine wate
    denbwnogas = (dbw_T_70 * np.exp((Ib_T_P_m - Ib_T_70_m)))

    # Step3: Calculate the vapor pressure of pure water

    T_array = T_array + 273

    # Coefficients Table 4_9
    a1Po = -7.85951783
    a2Po = 1.84408259
    a3Po = -11.7866497
    a4Po = 22.6807411
    a5Po = -15.9618719
    a6Po = 1.80122502

    # Ec 4.13
    vad = 1 - T_array / 674.096  ##Adimensional T EC.4.14
    # Ec 4.14
    Po = 22.064 * (np.exp((674.096 / T_array) * (
                a1Po * vad + a2Po * vad ** 1.5 + a3Po * vad ** 3 + a4Po * vad ** 3.5 + a5Po * vad ** 4 + a6Po * vad ** 7.5)))

    # Step4: Calculate the solubility of methane in pure water at (P,T)

    T_array = T_array - 273

    # Coefficients Table 4_10 EC 4.1
    a1A_T = 0
    a2A_T = -0.004462
    a3A_T = -0.06763
    a4A_T = 0
    a5A_T = 0

    A_T = (a1A_T * (T_array / 100) ** 2 + a2A_T * (T_array / 100) + a3A_T) / (
                a4A_T * (T_array / 100) ** 2 + a5A_T * (T_array / 100) + 1)

    a1B_T = -0.03602
    a2B_T = 0.18917
    a3B_T = 0.97242
    a4B_T = 0
    a5B_T = 0

    B_T = (a1B_T * (T_array / 100) ** 2 + a2B_T * (T_array / 100) + a3B_T) / (
                a4B_T * (T_array / 100) ** 2 + a5B_T * (T_array / 100) + 1)

    a1C_T = 0.6855
    a2C_T = -3.1992
    a3C_T = -3.7968
    a4C_T = 0.07711
    a5C_T = 0.2229

    C_T = (a1C_T * (T_array / 100) ** 2 + a2C_T * (T_array / 100) + a3C_T) / (
                a4C_T * (T_array / 100) ** 2 + a5C_T * (T_array / 100) + 1)

    # Ec 4.13 (Solubility of methane in pure water)
    mCH4pw = np.exp(A_T * (np.log(P_array - Po)) ** 2 + B_T * np.log(P_array - Po) + C_T)  ##Ec 4.15

    T_array = T_array + 273

    # Step5: Calculate the solubility of methane in brine water at (P,T)

    # Coefficients Table 4_11
    lambC1 = -0.80898
    lambC2 = 1.0827 * 10 ** (-3)
    lambC3 = 183.85
    lambC6 = 3.924 * 10 ** (-4)
    lambC10 = -1.97 * 10 ** (-6)

    # Ec 4.16
    lambCH4 = lambC1 + lambC2 * T_array + lambC3 / T_array + lambC6 * P_array + lambC10 * P_array ** 2  # EC4.16

    # Ec 4.18
    mCH4bw = mCH4pw * np.exp(-2 * lambCH4 * m - (-3.89 * 10 ** (-3)) * m ** 2)

    # Ec. 4.23 Specific Volume of methane-free
    Vbw = 1 / (denbwnogas)

    # Step6: Calculate the formation volume factor

    # Ec 4.36

    if unit == 1:
        Bw = ((1000 + m * 58.44) * Vbw + mCH4bw * VMCH4B) / ((1000 + m * 58.44) * 0.99681786)
        Bw = Bw * 1

    else:
        Bw = ((1000 + m * 58.44) * Vbw + mCH4bw * VMCH4B) / ((1000 + m * 58.44) * 0.99681786)

    return Bw


def comp_bw_nogas(P, T, salinity, unit=1):
    """This function calculate the compressibility of methane-free brine
    The user must select UNITS 1 to Field units and 2 to International System

    IF UNITS 1 is selected introduce:
        | Pressure (P) in PSI
        | Temperature (T) in ºF
        | Salinity in ppm
        | Compressibility returns in psi-1

    IF UNITS 2 is selected introduce:
        | Pressure (P) in MPa
        | Temperature (T) in ºC
        | Salinity in ppm
        | Compressibility returns in MPa-1"""

    if unit == 1:
        P_array = np.array(P)
        P_array = P_array * (6.89475729 * 10 ** -3.0)
        T_array = np.array(T)
        T_array = (T_array - 32) * 5 / 9
        m = np.array(salinity)
        m = 1000 * (m / 1000000) / (58.44 * (1 - (m / 1000000)))
    else:
        P_array = np.array(P)
        T_array = np.array(T)
        m = np.array(salinity)
        m = 1000 * (m / 1000000) / (58.44 * (1 - (m / 1000000)))

    # Step 1: Calculate the compressibility coefficients

    # Ec 4.1 and coefficients Table 4_6
    a1_EW = 4.221
    a2_EW = -3.478
    a3_EW = 6.221
    a4_EW = 0.5182
    a5_EW = -0.4405

    a_T_EW = (a1_EW * (T_array / 100) ** 2 + a2_EW * (T_array / 100) + a3_EW) / (
                a4_EW * (T_array / 100) ** 2 + a5_EW * (T_array / 100) + 1)

    a1_FW = -11.403
    a2_FW = 29.932
    a3_FW = 27.952
    a4_FW = 0.20684
    a5_FW = 0.3768

    a_T_FW = (a1_FW * (T_array / 100) ** 2 + a2_FW * (T_array / 100) + a3_FW) / (
            a4_FW * (T_array / 100) ** 2 + a5_FW * (T_array / 100) + 1)

    # Ec 4.1 and coefficients Table 4_8

    a1_EM = 0
    a2_EM = 0
    a3_EM = 0.1249
    a4_EM = 0
    a5_EM = 0

    a_T_EM = (a1_EM * (T_array / 100) ** 2 + a2_EM * (T_array / 100) + a3_EM) / (
            a4_EM * (T_array / 100) ** 2 + a5_EM * (T_array / 100) + 1)

    a1_FM3_2 = -0.617
    a2_FM3_2 = -0.747
    a3_FM3_2 = -0.4339
    a4_FM3_2 = 0
    a5_FM3_2 = 10.26

    a_T_FM3_2 = (a1_FM3_2 * (T_array / 100) ** 2 + a2_FM3_2 * (T_array / 100) + a3_FM3_2) / (
            a4_FM3_2 * (T_array / 100) ** 2 + a5_FM3_2 * (T_array / 100) + 1)

    a1_FM1 = 0
    a2_FM1 = 9.917
    a3_FM1 = 5.1128
    a4_FM1 = 0
    a5_FM1 = 3.892

    a_T_FM1 = (a1_FM1 * (T_array / 100) ** 2 + a2_FM1 * (T_array / 100) + a3_FM1) / (
                a4_FM1 * (T_array / 100) ** 2 + a5_FM1 * (T_array / 100) + 1)

    a1_FM1_2 = 0.0365
    a2_FM1_2 = -0.0369
    a3_FM1_2 = 0
    a4_FM1_2 = 0
    a5_FM1_2 = 0

    a_T_FM1_2 = (a1_FM1_2 * (T_array / 100) ** 2 + a2_FM1_2 * (T_array / 100) + a3_FM1_2) / (
                a4_FM1_2 * (T_array / 100) ** 2 + a5_FM1_2 * (T_array / 100) + 1)

    # Ec 4.7
    Eb_T_M = a_T_EW + a_T_EM * m
    # Ec 4.8
    FB_T_M = a_T_FW + a_T_FM3_2 * m ** (3 / 2) + a_T_FM1 * m + a_T_FM1_2 * m ** (1 / 2)

    # Step 2: Calculate the Compressibility of methane-free Brine at (P,T,m)
    # Ec 4.9  (Compressibility of methane-free Brine at (P,T,m))

    if unit == 1:
        temp = 1 / 70 * (1 / (Eb_T_M * (P_array / 70) + FB_T_M))
        temp = temp * 6.89475729 * 10 ** (-3.0)

    else:
        temp = 1 / 70 * (1 / (Eb_T_M * (P_array / 70) + FB_T_M))

    return temp
