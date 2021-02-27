# -*- coding: utf-8 -*-

"""
# Partie Bastien

deflectionAngleFromShock
machWaveAngle
shockAngleFromDeflection
maximumDeflectionAngle
shockWaveRange

# Partie Olivier

p_sur_pisentropique
T_sur_Tisentropique
rho_sur_rhoisentropique

# Partie Anthony

prandtlMeyerNuMuFromMach
prandtlMeyerMachFromNu
prandtlMeyerNuFromMach
prandltMeyerMuFromMach
maximumNu

# Partie Antoine

mach2frommach1
p2surp1frommach1
t2surt1frommach1
rho2surrho1frommach1
po2surpo1frommach1

# Tangent wedge method

p2surp1fromMachSh
rho2surrho1fromMachSh
t2surt1fromMachSh

# Modified Newtonian Technique

p2surp1fromM0De
cpFromM0p2surp1

# US standard atmosphere extended to 1000 km

geopotentialFromGeometric   km -> km'
geometricFromGeopotential   km' -> km
temperatureFromAlt          km -> K
getGFromAlt                 km -> m/s²
pressureFromAlt             km -> Pa
densityFromAlt              km -> kg/m^3
soundSpeedFromAlt           km -> m/s
"""

from scipy.interpolate import interp1d
from scipy.optimize import fminbound, bisect
import numpy as np

gamma = 1.4
R_T = 6356.766  # km
g_0 = 9.80665  # m/s²
k_B = 1.380622e-23  # J/K
R = 287.053  # J/(kg.K)
Rstar = 0.0083144621  # kg(km/s)²/(kmol.K)
M_0 = 28.9644  # kg/kmol

M_H = 1.00794  # kg/kmol
M_He = 4.002602  # kg/kmol
M_O = 15.9994  # kg/kmol

P_0 = 101325.0  # Pa
P_1 = 22632.06  # Pa
P_2 = 5474.889  # Pa
P_3 = 868.0187  # Pa
P_4 = 110.9063  # Pa
P_5 = 66.93887  # Pa
P_6 = 3.956420  # Pa

L_M0 = -6.5  # K/km'
L_M2 = +1.0  # K/km'
L_M3 = +2.8  # K/km'
L_M5 = -2.8  # K/km'
L_M6 = -2.0  # K/km'
L_K9 = 12  # K/km

T_M0 = 288.15  # K
T_M1 = 216.65  # K
T_M3 = 228.65  # K
T_M4 = 270.65  # K
T_M6 = 214.65  # K
T_7 = 186.8673  # K
T_c = 263.1905  # K
T_9 = 240  # K
T_10 = 360.0  # K
T_inf = 1000  # K

r_H1000 = 0.0395115262
r_He1000 = 0.930515916
r_O1000 = 0.0299725578

minM, maxM = 1, 100
table = {}

fileDSM = "DSMdata.csv"  # Deflection-Shock-Mach (θ-β-M)
fileIFP = "IFPdata.csv"  # Isentropic Flow Properties
fileMNuMu = "MNuMudata.csv"  # Prandtl-Meyer
fileNSP = "Shockdroitdata.csv"  # Normal Shock Properties


def writeTable(file, data):
    np.savetxt(file, data, delimiter=",")


def readTable(file):
    return np.genfromtxt(file, delimiter=",")


def deflectionAngleFromShock(delta, M):
    return np.arctan2(2 * np.cos(delta) * ((M * np.sin(delta))**2 - 1),
                      np.sin(delta) * (M**2 * (gamma + np.cos(2 * delta)) + 2))


def machWaveAngle(M):
    return np.arcsin(1 / M)


def shockAngleFromDeflection(theta, M):
    # Find angle using bisection
    assert theta >= 0

    if theta <= maximumDeflectionAngle(M):
        return bisect(lambda d: (deflectionAngleFromShock(d, M) - theta),
                      *shockWaveRange(M))
    raise "Deflection angle too big. No attached shock."


def maximumDeflectionAngle(M):
    return table["MDA"](M)


def shockWaveRange(M):
    return table["SWAmin"](M), table["SWAmax"](M)

# Partie Olivier


def p_sur_pisentropique(M):
    return (1 + ((gamma - 1) * M**2 / 2))**(-gamma / (gamma - 1))


def T_sur_Tisentropique(M):
    return (1 + ((gamma - 1) * M**2 / 2))**(-1)


def rho_sur_rhoisentropique(M):
    return (1 + ((gamma - 1) * M**2 / 2))**(-1 / (gamma - 1))

# Partie Anthony


def prandtlMeyerNuMuFromMach(M):
    return (np.sqrt((gamma + 1) / (gamma - 1)) * np.arctan(
            np.sqrt((gamma - 1) / (gamma + 1) * (M**2 - 1))) - np.arctan(
            np.sqrt(M**2 - 1))), np.arcsin(1 / M)


@np.vectorize  # Decorator (machFromNu = np.vectorize(machFromNu))
def prandtlMeyerMachFromNu(Nu):
    return bisect(lambda m: (table["Nu"](m) - Nu), minM, maxM)


def prandtlMeyerNuFromMach(M):
    return table["Nu"](M)


def prandltMeyerMuFromMach(M):
    return np.arcsin(1 / M)  # table["Mu"](M)


def maximumNu():
    # Limited by Mach table not theory
    # approaches real limitation as maxM increases to inf:
    # nu_max = np.pi * (np.sqrt((gamma + 1) / (gamma - 1)) - 1) / 2
    return table["maxNu"]

# Partie Antoine


def mach2frommach1(M):  # calcul of mach2 from mach1
    return (np.sqrt((M**2 * (gamma - 1) + 2)
                    / (M**2 * 2 * gamma - (gamma - 1))))


def p2surp1frommach1(M):  # calcul of p2 with p1 and mach1
    return ((2 * gamma * M**2) / (gamma + 1)) - ((gamma - 1) / (gamma + 1))


def t2surt1frommach1(M):  # Calcul of t2/t1
    return (((1 + (gamma - 1) / 2 * M**2)
             * (2 * gamma * M**2 / (gamma - 1) - 1))
            / (M**2 * ((2 * gamma) / (gamma - 1) + (gamma - 1) / 2)))


def rho2surrho1frommach1(M):  # Calcul de rho2/rho1
    return ((gamma + 1) * M**2) / ((gamma - 1) * M**2 + 2)


def po2surpo1frommach1(M):
    return ((((gamma + 1) * M**2 / 2)
             / (1 + (gamma - 1) / 2 * M**2))**(gamma / (gamma - 1))
            * (1 / ((2 * gamma) * M**2 / (gamma + 1)
                    - (gamma - 1) / (gamma + 1)))**(1 / (gamma - 1)))

# Tangent wedge method


def p2surp1fromMachSh(M, sAngle):  # Calcul de p2/p1
    return 1 + 2 * gamma * (M**2 * np.sin(sAngle)**2 - 1) / (gamma + 1)


def rho2surrho1fromMachSh(M, sAngle):  # Calcul de rho2/rho1
    return (gamma + 1) * M**2 * np.sin(sAngle)**2 / \
        ((gamma - 1) * M**2 * np.sin(sAngle)**2 + 2)


def t2surt1fromMachSh(M, sAngle):  # Calcul de t2/t1
    return p2surp1fromMachSh(M, sAngle) / rho2surrho1fromMachSh(M, sAngle)


# Modified Newtonian Technique

def p2surp1fromM0De(M0, dAngle):
    return 1 + gamma * M0**2 * np.sin(dAngle)**2


def cpFromM0p2surp1(M0, p2surp1):
    return 2 * (p2surp1 - 1) / (gamma * M0**2)


def interp(x, y, kind='cubic'):
    return interp1d(x, y, kind)


def main(plt=False):
    """ Generate tables """

    n = 200
    minMlog = np.log10(minM)
    maxMlog = np.log10(maxM)
    moyMlog = (minMlog + maxMlog) / 2
    listM1 = np.logspace(minMlog, maxMlog, n)
    listMDA = np.empty(n)  # Maximum Deflection Angle
    listSWAmin = np.empty(n)  # Shock Wave Angle at minimum
    listSWAmax = np.empty(n)  # Shock Wave Angle at maximum
    listMu = np.empty(n)  # Mu Prandlt Meyer
    listNu = np.empty(n)  # Nu Prandlt Meyer
    listP2 = p2surp1frommach1(listM1)  # Shock Droit Pression
    listT2 = t2surt1frommach1(listM1)  # Shock Droit Température
    listSM2 = mach2frommach1(listM1)  # Shock Droit Mach2
    listpo2 = po2surpo1frommach1(listM1)  # Shock Droit Pression02
    listrho2 = rho2surrho1frommach1(listM1)  # Shock Droit rho2

    iterM = enumerate(listM1)
    next(iterM)  # Skip first iteration
    listMDA[0] = 0
    listSWAmin[0] = np.pi / 2
    listSWAmax[0] = np.pi / 2
    listNu[0] = 0
    listMu[0] = np.pi / 2
    for i, M in iterM:
        xopt, fval, ierr, numfunc = fminbound(
            lambda t: -deflectionAngleFromShock(t, M), 0, np.pi / 2,
            full_output=True)
        assert ierr == 0
        listMDA[i] = -fval
        listSWAmax[i] = xopt

        x0, r = bisect(deflectionAngleFromShock, 0, xopt, args=M,
                       full_output=True)
        assert r.converged
        listSWAmin[i] = x0
        listNu[i], listMu[i] = prandtlMeyerNuMuFromMach(M)

    # Partie Olivier
    n1 = 200
    n2 = 100
    listM21 = np.logspace(minMlog, moyMlog, n1, endpoint=False)  # Mach
    listM22 = np.logspace(moyMlog, maxMlog, n2)
    listM2 = np.concatenate((listM21, listM22))

    # pressure/isentropic pressure
    listp_sur_pi = p_sur_pisentropique(listM2)
    # temperature/isentropic temperature
    listT_sur_Ti = T_sur_Tisentropique(listM2)
    # density/isentropic density
    listrho_sur_rhoi = rho_sur_rhoisentropique(listM2)

    writeTable(fileDSM, (listM1, listMDA, listSWAmin, listSWAmax))
    table["MDA"] = interp(listM1, listMDA)
    table["SWAmin"] = interp(listM1, listSWAmin)
    table["SWAmax"] = interp(listM1, listSWAmax)

    writeTable(fileIFP, (listM2, listp_sur_pi,
                         listT_sur_Ti, listrho_sur_rhoi))
    table["p_sur_pi"] = interp(listM2, listp_sur_pi)
    table["T_sur_Ti"] = interp(listM2, listT_sur_Ti)
    table["rho_sur_rhoi"] = interp(listM2, listrho_sur_rhoi)

    writeTable(fileMNuMu, (listM1, listNu, listMu))
    table["Nu"] = interp(listM1, listNu)
    table["Mu"] = interp(listM1, listMu)
    table["maxNu"] = max(listNu)

    writeTable(fileNSP, (listM1, listP2, listSM2, listT2, listrho2, listpo2))
    table["P2"] = interp(listM1, listP2)
    table["SM2"] = interp(listM1, listSM2)
    table["T2"] = interp(listM1, listT2)
    table["rho2"] = interp(listM1, listrho2)
    table["po2"] = interp(listM1, listpo2)

    # Tables Prandtl-Meyer

    if plt:
        from matplotlib.ticker import ScalarFormatter

        fig1, ax1 = plt.subplots()
        ax1.plot(listM1, listMDA * 180 / np.pi,
                 label="Maximum deflection")
        ax1.plot(listM1, listSWAmin * 180 / np.pi,
                 label="Shock wave minimum")
        ax1.plot(listM1, listSWAmax * 180 / np.pi,
                 label="Shock wave maximum")
        ax1.legend()
        ax1.set_title("Shock wave angle at maximum and at no deflection"
                      " (γ = " + str(gamma) + ")")
        ax1.set_xlabel("Mach")
        ax1.set_ylabel("Angle (°)")
        ax1.set_xscale("log")
        ax1.grid(True, which="both")
        ax1.set_xticks([1, 2, 3, 5, 10, 20, 30, 50, 100])
        yticks = list(range(0, 90, 20))
        yticks.append(round(listSWAmax[-1] * 180 / np.pi))
        yticks.append(round(listMDA[-1] * 180 / np.pi))
        ax1.set_yticks(yticks)
        ax1.xaxis.set_major_formatter(ScalarFormatter())

        fig2, ax2 = plt.subplots()
        ax2.plot(listM2, listp_sur_pi, label="P/P\u2080")
        ax2.plot(listM2, listT_sur_Ti, label="T/P\u2080")
        ax2.plot(listM2, listrho_sur_rhoi, label="ρ/ρ\u2080")
        ax2.legend()
        ax2.set_title("Isentropic flow properties (γ = " + str(gamma) + ")")
        ax2.set_xlabel("Mach")
        ax2.set_xscale("log")
        ax2.grid(True, which="both")
        ax2.set_xticks([1, 2, 3, 5, 10, 20, 30, 50, 100])
        ax2.xaxis.set_major_formatter(ScalarFormatter())

        fig3, ax3 = plt.subplots()
        ax3.plot(listM1, listMu * 180 / np.pi, label="Mu")
        ax3.plot(listM1, listNu * 180 / np.pi, label="Nu")
        ax3.legend()
        ax3.set_title("Prandtl-Meyer function and Mach angle"
                      " (γ = " + str(gamma) + ")")
        ax3.set_xlabel("Mach")
        ax3.set_ylabel("Angle (°)")
        ax3.set_xscale("log")
        ax3.grid(True, which="both")
        ax3.set_xticks([1, 2, 3, 5, 10, 20, 30, 50, 100])
        ax3.xaxis.set_major_formatter(ScalarFormatter())

        fig4, ax4 = plt.subplots()
        ax4.plot(listM1, listP2, label="p\u2082/p\u2081")
        ax4.plot(listM1, listT2, label="T\u2082/T\u2081")
        ax4.legend()
        ax4.set_title("Normal Shock Properties"
                      " (γ = " + str(gamma) + ")")
        ax4.set_xlabel("Mach")
        ax4.set_ylabel("Ratio")
        ax4.set_xscale("log")
        ax4.set_yscale("log")
        ax4.grid(True, which="both")
        ax4.set_xticks([1, 2, 3, 5, 10, 20, 30, 50, 100])
        ax4.xaxis.set_major_formatter(ScalarFormatter())

        fig5, ax5 = plt.subplots()
        ax5.plot(listM1, listSM2, label="M\u2082")
        ax5.plot(listM1, listpo2, label="po\u2082/po\u2081")
        ax5.plot(listM1, listrho2, label="rho\u2082/rho\u2081")
        ax5.legend()
        ax5.set_title("Normal Shock Properties"
                      " (γ = " + str(gamma) + ")")
        ax5.set_xlabel("Mach")
        ax5.set_ylabel("Mach")
        ax5.set_xscale("log")
        ax5.grid(True, which="both")
        ax5.set_xticks([1, 2, 3, 5, 10, 20, 30, 50, 100])
        ax5.xaxis.set_major_formatter(ScalarFormatter())

        plt.show()


coefPressure = {
    86: [0.000000, 2.159582E-06, -4.836957E-04, -0.1425192, 13.47530],
    91: [0.000000, 3.304895E-05, -0.009062730, 0.6516698, -11.03037],
    100: [0.000000, 6.693926E-05, -0.01945388, 1.719080, -47.75030],
    110: [0.000000, -6.539316E-05, 0.02485568, -3.223620, 135.9355],
    120: [2.283506E-07, -1.343221E-04, 0.02999016, -3.055446, 113.5764],
    150: [1.209434E-08, -9.692458E-06, 0.003002041, -0.4523015, 19.19151],
    200: [8.113942E-10, -9.822568E-07, 4.687616E-04, -0.1231710, 3.067409],
    300: [9.814674E-11, -1.654439E-07, 1.148115E-04, -0.05431334, -2.011365],
    500: [-7.835161E-11, 1.964589E-07, -1.657213E-04, 0.04305869, -14.77132],
    750: [2.813255E-11, -1.120689E-07, 1.695568E-04, -0.1188941, 14.56718]
}

coefDensity = {
    86: [0.000000, -3.322622E-06, 9.111460E-04, -0.2609971, 5.944694],
    91: [0.000000, 2.873405E-05, -0.008492037, 0.6541179, -23.62010],
    100: [-1.240774E-05, 0.005162063, -0.8048342, 55.55996, -1443.338],
    110: [0.00000, -8.854164E-05, 0.03373254, -4.390837, 176.5294],
    120: [3.661771E-07, -2.154344E-04, 0.04809214, -4.884744, 172.3597],
    150: [1.906032E-08, -1.527799E-05, 0.004724294, -0.6992340, 20.50921],
    200: [1.199282E-09, -1.451051E-06, 6.910474E-04, -0.1736220, -5.321644],
    300: [1.140564E-10, -2.130756E-07, 1.570762E-04, -0.07029296, -12.89844],
    500: [8.105631E-12, -2.358417E-09, -2.635110E-06, -0.01562608, -20.02246],
    750: [-3.701195E-12, -8.608611E-09, 5.118829E-05, -0.06600998, -6.137674]
}


def geopotentialFromGeometric(z):  # km -> km'
    return R_T * z / (R_T + z)


def geometricFromGeopotential(h):  # km' -> km
    return R_T * h / (R_T - h)


@np.vectorize
def temperatureFromAlt(z):  # km -> K
    if z > 120:  # Exosphere / Thermosphere
        xi = (z - 120) * (R_T + 120) / (R_T + z)
        return T_inf - (T_inf - T_10) * np.exp(-0.01875 * xi)
    elif z > 110:
        return T_9 + L_K9 * (z - 110)
    elif z > 91:
        return T_c - 76.3232 * np.sqrt(1 - ((91 - z) / 19.9429)**2)
    elif z > 86:  # Mesosphere
        return T_7

    h = geopotentialFromGeometric(z)
    if z > 80:
        return (T_M6 + L_M6 * (h - 71)) * (5028080 - 351 * z) / 5000000
    elif h > 71:
        return T_M6 + L_M6 * (h - 71)
    elif h > 51:
        return T_M4 + L_M5 * (h - 51)
    elif h > 47:
        return T_M4
    elif h > 32:  # Stratosphere
        return T_M3 + L_M3 * (h - 32)
    elif h > 20:
        return T_M1 + L_M2 * (h - 20)
    elif h > 11:
        return T_M1
    else:  # Troposphere
        return T_M0 + L_M0 * h


def getGFromAlt(z):  # km -> m/s²
    return g_0 * (R_T / (R_T + z))**2


@np.vectorize
def pressureFromAlt(z):  # km -> Pa
    # Not accurate after 1000km (Guess from molecular weight)

    if z > 1000:  # Exosphere / Thermosphere
        i = list(coefPressure.keys())[-1]
        maxAltP = np.exp(coefPressure[i][0] * 1000**4
                         + coefPressure[i][1] * 1000**3
                         + coefPressure[i][2] * 1000**2
                         + coefPressure[i][3] * 1000
                         + coefPressure[i][4])  # Pa
        n_tot = maxAltP / (k_B * T_inf)  # m^-3

        g = getGFromAlt(z) / 1000  # km/s²
        gRsT = g / (Rstar * T_inf)  # kmol/(kg.km)

        n_H = n_tot * r_H1000 * np.exp((1000 - z) * M_H * gRsT)  # m^-3
        n_He = n_tot * r_He1000 * np.exp((1000 - z) * M_He * gRsT)  # m^-3
        n_O = n_tot * r_O1000 * np.exp((1000 - z) * M_O * gRsT)  # m^-3

        return (n_H + n_He + n_O) * k_B * T_inf
    elif z > 86:
        k = list(coefPressure.keys())
        i, n = 1, len(k)
        while i < n and k[i] < z:
            i += 1
        i = k[i - 1]
        return np.exp(coefPressure[i][0] * z**4 + coefPressure[i][1] * z**3
                      + coefPressure[i][2] * z**2 + coefPressure[i][3] * z
                      + coefPressure[i][4])

    h = geopotentialFromGeometric(z)
    GMR = 34.1632  # km'/K
    if h > 71:  # Mesosphere
        return P_6 * (T_M6 / (T_M6 + L_M6 * (h - 71)))**(GMR / L_M6)
    elif h > 51:
        return P_5 * (T_M4 / (T_M4 + L_M5 * (h - 51)))**(GMR / L_M5)
    elif h > 47:
        return P_4 * np.exp(-GMR * (h - 47) / T_M4)
    elif h > 32:  # Stratosphere
        return P_3 * (T_M3 / (T_M3 + L_M3 * (h - 32)))**(GMR / L_M3)
    elif h > 20:
        return P_2 * (T_M1 / (T_M1 + (h - 20)))**GMR
    elif h > 11:
        return P_1 * np.exp(-GMR * (h - 11) / T_M1)
    else:  # Troposphere
        return P_0 * (T_M0 / (T_M0 + L_M0 * h))**(GMR / L_M0)


@np.vectorize
def densityFromAlt(z):  # km -> kg/m^3
    # Not accurate after 1000km (Polynomial extrapolation)
    if z < 80:
        return pressureFromAlt(z) / (R * temperatureFromAlt(z))
    elif z < 86:
        h = geopotentialFromGeometric(z)
        return pressureFromAlt(z) / (R * (T_M6 + L_M6 * (h - 71)))

    k = list(coefPressure.keys())
    i, n = 1, len(k)
    while i < n and k[i] < z:
        i += 1
    i = k[i - 1]
    return np.exp(coefDensity[i][0] * z**4 + coefDensity[i][1] * z**3
                  + coefDensity[i][2] * z**2 + coefDensity[i][3] * z
                  + coefDensity[i][4])


def soundSpeedFromTemp(t):  # K -> m/s
    return 1000 * np.sqrt(gamma * Rstar * t / M_0)


def soundSpeedFromAlt(z):  # km -> m/s
    # Not applicable for high altitudes (> 86km)
    # (temperatureFromAlt returns T_K, needs T_M)
    return soundSpeedFromTemp(temperatureFromAlt(z))


def plotAtmo(plt):
    # Figure_P

    z = np.linspace(0, 100, 500)
    fig1, ax1 = plt.subplots()

    ax1.plot(pressureFromAlt(z), z)
    ax1.set_title("Pression en fonction de l'altitude")
    ax1.set_xlabel("Pression (Pa)")
    ax1.set_ylabel("Altitude (km)")
    ax1.grid(True, which="both")

    # Figure_T

    z = np.linspace(0, 100, 500)
    fig2, ax2 = plt.subplots()

    ax2.plot(temperatureFromAlt(z), z)
    ax2.set_title("Température en fonction de l'altitude")
    ax2.set_xlabel("Température (K)")
    ax2.set_ylabel("Altitude (km)")
    ax2.grid(True, which="both")

    # Figure_D

    z = np.linspace(0, 100, 500)
    fig3, ax3 = plt.subplots()

    ax3.plot(densityFromAlt(z), z)
    ax3.set_title("Densité en fonction de l'altitude")
    ax3.set_xlabel("Densité (kg/m\u00B3)")
    ax3.set_ylabel("Altitude (km)")
    ax3.set_xscale("log")
    ax3.grid(True, which="both")

    # Figure_S

    z = np.linspace(0, 100, 500)
    fig4, ax4 = plt.subplots()

    ax4.plot(soundSpeedFromAlt(z), z)
    ax4.set_title("Vitesse du son en fonction de l'altitude")
    ax4.set_xlabel("Vitesse (m/s)")
    ax4.set_ylabel("Altitude (km)")
    ax4.grid(True, which="both")

    plt.show()


if __name__ == "__main__":
    import matplotlib.pyplot as plt
    main(plt)

    # mach = np.linspace(3, 25, 200)
    # fig1, ax1 = plt.subplots(figsize=(10, 7))

    # p1 = pressureFromAlt(10)
    # t1 = temperatureFromAlt(10)
    # rho1 = densityFromAlt(10)
    # c = soundSpeedFromAlt(10)

    # theta = 5  # 5, 15, 30
    # dAngle = theta * np.pi / 180
    # dataC = []
    # dataTPP = []
    # dataDT = []
    # dataCD = []
    # dataNS = []
    # dataNM = []

    # for M1 in mach:
    #     # C
    #     sAngleC = shockAngleFromDeflection(dAngle, M1)

    #     Mn1 = M1 * np.sin(sAngleC)
    #     Mn2 = mach2frommach1(Mn1)
    #     M2C = Mn2 / np.sin(sAngleC - dAngle)

    #     Mt1 = M1 * np.cos(sAngleC)
    #     Vt1 = c * Mt1

    #     t2C = M_0 * Vt1**2 / (1e6 * gamma * Rstar * (M2C**2 - Mn2**2))

    #     Vn1 = c * Mn1
    #     Vn2 = c * Mn2
    #     rho2C = rho1 * Vn1 / Vn2

    #     p2C = p1 + rho1 * Vn1**2 - rho2C * Vn2**2
    #     dataC.append(sAngleC)

    #     # TPP
    #     sAngleTTP = dAngle * ((gamma + 1) / 4
    #                           + np.sqrt((gamma + 1)**2 / 16
    #                                     + 1 / (M1**2 * dAngle**2)))

    #     p2TTP = p1 + p1 * 2 * gamma * \
    #         (M1**2 * np.sin(sAngleTTP)**2 - 1) / (gamma + 1)

    #     p2TTP2 = 2 * gamma * M1**2 * np.sin(sAngleTTP)**2 / (gamma + 1) \
    #         - (gamma - 1) / (gamma + 1)

    #     MSA2 = M1**2 * sAngleTTP**2
    #     M2TTP = np.sqrt((gamma + 1)**2 * M1**2 * MSA2
    #                     - 4 * (MSA2 - 1) * (gamma * MSA2 + 1)
    #                     / ((2 * gamma * MSA2 - gamma + 1)
    #                         * ((gamma - 1) * MSA2 + 2)))

    #     rho2TTP = rho1 * (gamma + 1) * M1**2 * np.sin(sAngleTTP)**2 \
    #         / ((gamma - 1) * M1**2 * np.sin(sAngleTTP)**2 + 2)

    #     t2TTP = t1 * (p2TTP / p1) / (rho2TTP / rho1)
    #     t2TTP2 = t1 * (p2TTP2 / p1) / (rho2TTP / rho1)

    #     dataTPP.append(sAngleTTP)

    #     # DT
    #     sAngleDT = shockAngleFromDeflection(dAngle, M1)

    #     Mn1 = M1 * np.sin(sAngleDT)
    #     Mn2 = mach2frommach1(Mn1)
    #     M2DT = Mn2 / np.sin(sAngleDT - dAngle)

    #     p2DT = p1 + p1 * 2 * gamma * \
    #         (M1**2 * np.sin(sAngleDT)**2 - 1) / (gamma + 1)

    #     rho2DT = rho1 * rho2surrho1frommach1(Mn1)

    #     t2DT = t1 * p2DT * rho1 / (p1 * rho2DT)

    #     dataDT.append(sAngleDT)

    #     # CD
    #     sAngleCD = shockAngleFromDeflection(dAngle, M1)

    #     Mn1 = M1 * np.sin(sAngleCD)
    #     Mn2 = mach2frommach1(Mn1)
    #     M2CD = Mn2 / np.sin(sAngleCD - dAngle)

    #     nu2 = prandtlMeyerNuFromMach(M2CD)
    #     delta = dAngle - dAngle  # θ_i - θ_c
    #     MiCD = prandtlMeyerMachFromNu(delta + nu2)

    #     p2CD = p1 * p2surp1frommach1(Mn1)
    #     # piCD = p2CD * p_sur_pisentropique(MiCD)
    #     t2CD = t1 * t2surt1frommach1(Mn1)
    #     rho2CD = rho1 * rho2surrho1frommach1(Mn1)

    #     dataCD.append(sAngleCD)

    #     # NS
    #     p2NS = p1 + p1 * gamma * M1**2 * np.sin(dAngle)**2
    #     Vn = M1 * c * np.sin(dAngle)
    #     M2NS = Vn / c  # No temperature change?
    #     dataNS.append(M2NS)

    #     NM
    #     p2NM = p1 + rho1 * (M1 * c)**2 * np.sin(dAngle)**2
    #     Vn = M1 * c * np.sin(dAngle)
    #     M2NM = Vn / c  # No temperature change?
    #     dataNM.append(M2NM)

    # dataC = np.array(dataC)  # * 180 / np.pi
    # dataTPP = np.array(dataTPP)  # * 180 / np.pi
    # dataDT = np.array(dataDT)  # * 180 / np.pi
    # dataCD = np.array(dataCD)  # * 180 / np.pi
    # dataNS = np.array(dataNS)
    # dataNM = np.array(dataNM)

    # n = 2
    # ax1.plot(mach, mach, "grey", label=r"$\mathrm{M_\infty}$",
    #          linewidth=n)
    # ax1.plot(mach, dataC, label="Classique", linewidth=1.5 * n)
    # ax1.plot(mach, dataTPP, label="TPP", linewidth=n)
    # ax1.plot(mach, dataDT, label="Dièdre tangent", linewidth=1.5 * n,
    #          linestyle="dashed")
    # ax1.plot(mach, dataCD, label="Choc détente", linewidth=2 * n,
    #          linestyle="dotted")
    # ax1.plot(mach, dataNS, label="Newton simplifié", linewidth=n)
    # ax1.plot(mach, dataNM, label="Newton modifié", linewidth=2 * n,
    #          linestyle="dotted")

    # ax1.legend()
    # ax1.set_title("Angle du choc pour différentes méthodes"
    #               " (10km d'altitude, θ=%i°)" % theta)
    # ax1.set_xlabel("Mach")
    # ax1.set_ylabel("Angle (°)")
    # ax1.grid(True, which="both")

    # plt.show()
else:
    try:
        listM, listMDA, listSWAmin, listSWAmax = readTable(fileDSM)
        table["MDA"] = interp(listM, listMDA)
        table["SWAmin"] = interp(listM, listSWAmax)
        table["SWAmax"] = interp(listM, listSWAmin)

        listM, listPsurPi, listTsurTi, listRhoSurRhoi = readTable(fileIFP)
        table["p_sur_pi"] = interp(listM, listPsurPi)
        table["T_sur_Ti"] = interp(listM, listTsurTi)
        table["rho_sur_rhoi"] = interp(listM, listRhoSurRhoi)

        listM, listNu, listMu = readTable(fileMNuMu)
        table["Nu"] = interp(listM, listNu)
        table["Mu"] = interp(listM, listMu)
        table["maxNu"] = max(listNu)

        listM, listP2, listSM2, listT2, listrho2, listpo2 = readTable(fileNSP)
        table["P2"] = interp(listM, listP2)
        table["SM2"] = interp(listM, listSM2)
        table["T2"] = interp(listM, listT2)
        table["rho2"] = interp(listM, listrho2)
        table["po2"] = interp(listM, listpo2)
    except OSError as e:
        print(e)
        main()
